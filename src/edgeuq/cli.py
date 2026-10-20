"""Command line: `edgeuq shift`, `edgeuq jitter`, `edgeuq rescan`, `edgeuq demo`."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from . import __version__
from .jitter import run_jitter, sigma_from_median
from .constants import N_SAMPLES, RADIAL_MM, SEED0
from .io import extract_endpoints, load_nodes, save_matrix
from .rescan import validate
from .shift import compute_u


def _json(path: Path, rec: dict) -> None:
    path.write_text(json.dumps(rec, indent=2, default=float) + "\n")
    print(json.dumps(rec, indent=2, default=float))


def _inputs(args) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    args.out.mkdir(parents=True, exist_ok=True)
    nodes, aff = load_nodes(args.nodes)
    n_nodes = args.n_nodes or int(nodes.max())
    cache = args.out / "endpoints.npy"
    if cache.is_file():
        ends = np.load(cache)
    else:
        ends = extract_endpoints(args.tck)
        np.save(cache, ends)
    return ends.astype(np.float32), nodes, aff, n_nodes


def _provenance(args, **extra) -> dict:
    return {
        "edgeuq_version": __version__,
        "tck": str(args.tck),
        "nodes": str(args.nodes),
        **extra,
    }


def cmd_shift(args) -> int:
    ends, nodes, aff, n_nodes = _inputs(args)
    shifts = args.out / "shifts"
    shifts.mkdir(exist_ok=True)
    res = compute_u(ends, nodes, aff, n_nodes, args.radius, cache_dir=shifts)
    save_matrix(args.out / "C0.csv", res["C0"], fmt="%d")
    save_matrix(args.out / "U.csv", res["U"])
    save_matrix(args.out / "Ucal.csv", res["Ucal"])
    save_matrix(args.out / "e_1mm.csv", res["e_1mm"])
    save_matrix(args.out / "mask.csv", res["mask"].astype(int), fmt="%d")
    _json(args.out / "shift.json", _provenance(args) | res["summary"])
    return 0


def cmd_jitter(args) -> int:
    if (args.sigma is None) == (args.leftover_median is None):
        raise SystemExit("give exactly one of --sigma or --leftover-median")
    sigma = args.sigma if args.sigma is not None else sigma_from_median(args.leftover_median)
    ends, nodes, aff, n_nodes = _inputs(args)

    def progress(s: int, n: int) -> None:
        print(f"[jitter] sample {s}/{n}", file=sys.stderr, flush=True)

    res = run_jitter(ends, nodes, aff, n_nodes, sigma, args.samples, args.seed, args.radius, progress)
    save_matrix(args.out / "F.csv", res["F"])
    save_matrix(args.out / "Esoft.csv", res["Esoft"])
    save_matrix(args.out / "P_ge_floor.csv", res["P_ge_floor"])
    if args.save_samples:
        np.save(args.out / "Cs.npy", res["Cs"])
    extra = {"leftover_median_mm": args.leftover_median} if args.leftover_median is not None else {}
    _json(args.out / "jitter.json", _provenance(args, **extra) | res["summary"])
    return 0


def cmd_rescan(args) -> int:
    pairs = [(Path(r), Path(a), Path(b)) for r, a, b in args.pair]
    rec = {"edgeuq_version": __version__} | validate(pairs)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    _json(args.out, rec)
    return 0


def cmd_demo(args) -> int:
    from .demo import write_demo

    p = write_demo(args.out / "inputs", n=args.streamlines)
    run = args.out / "run-1"
    common = ["--tck", str(p["tck"]), "--nodes", str(p["nodes"]), "--out", str(run)]
    main(["shift", *common])
    main(["jitter", *common, "--sigma", "0.5", "--samples", "8"])
    return main(["rescan", "--pair", str(run), str(p["count_a"]), str(p["count_b"]), "--out", str(args.out / "rescan.json")])


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="edgeuq", description="Edge-level uncertainty for streamline-count connectomes.")
    ap.add_argument("--version", action="version", version=f"edgeuq {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True, metavar="{shift,jitter,rescan,demo}")

    def add_inputs(p: argparse.ArgumentParser) -> None:
        p.add_argument("--tck", type=Path, required=True, help="tractogram (.tck or .tck.gz), scanner mm")
        p.add_argument("--nodes", type=Path, required=True, help="node image in the same space (.nii/.nii.gz, or .mif with mrconvert)")
        p.add_argument("--out", type=Path, required=True, help="output directory")
        p.add_argument("--n-nodes", type=int, default=None, help="number of nodes (default: max label; set it, e.g. 78 for DKT, if a label may be absent)")
        p.add_argument("--radius", type=float, default=RADIAL_MM, help=f"radial search mm (protocol: {RADIAL_MM:g})")

    p = sub.add_parser("shift", aliases=["u"], help="shift sensitivity: U and calibrated U from planted endpoint shifts")
    add_inputs(p)
    p.set_defaults(func=cmd_shift)

    p = sub.add_parser("jitter", aliases=["c2"], help="jitter propagation: leakage F, soft connectome, P(count >= floor)")
    add_inputs(p)
    p.add_argument("--sigma", type=float, help="per-axis endpoint jitter SD in mm")
    p.add_argument("--leftover-median", type=float, help="measured median leftover displacement in mm (sigma = median / 1.5382)")
    p.add_argument("--samples", type=int, default=N_SAMPLES)
    p.add_argument("--seed", type=int, default=SEED0)
    p.add_argument("--save-samples", action="store_true", help="also write per-sample counts Cs.npy")
    p.set_defaults(func=cmd_jitter)

    p = sub.add_parser("rescan", aliases=["retest"], help="rescan validation: score run-1 U/Ucal/F against run-1 vs run-2 change")
    p.add_argument("--pair", nargs=3, action="append", required=True, metavar=("RUN1_DIR", "COUNT_RUN1", "COUNT_RUN2"))
    p.add_argument("--out", type=Path, required=True, help="output JSON")
    p.set_defaults(func=cmd_rescan)

    p = sub.add_parser("demo", help="run shift, jitter and rescan on synthetic data (smoke test)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--streamlines", type=int, default=40_000)
    p.set_defaults(func=cmd_demo)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
