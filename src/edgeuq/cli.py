"""Command line: `edgeuq shift`, `edgeuq jitter`, `edgeuq rescan`, `edgeuq instability`, `edgeuq demo`."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from . import __version__
from .jitter import run_jitter, sigma_from_median
from .constants import COUNT_FLOOR, N_SAMPLES, RADIAL_MM, RELCHANGE_CUT, SEED0
from .instability import evaluate_leftover, score_hybrid, score_neg_log_count, score_one_minus_p10
from .io import extract_endpoints, load_matrix, load_nodes, save_matrix
from .leftover_cal import (
    StratumMismatchError,
    apply_leftover_map,
    attach_stratum_metadata,
    evaluate_PL,
    fit_PL,
    leftover_y,
    load_leftover_map,
    loso_PL,
)
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


def cmd_instability(args) -> int:
    """Ses-1 leftover-instability scores (Dir 5). Separate from shift/jitter; not hat{U}."""
    args.out.mkdir(parents=True, exist_ok=True)
    c0 = load_matrix(args.count_run1)
    L = score_neg_log_count(c0)
    save_matrix(args.out / "L.csv", L)
    save_matrix(args.out / "L_neg_log1p_C0.csv", L)
    summary: dict = {
        "edgeuq_version": __version__,
        "count_run1": str(args.count_run1),
        "primary": "L_neg_log1p_C0",
        "not_hatU": True,
        "distinct_from_shift_jitter": True,
        "formula": "L = -log1p(C0)",
    }
    p10 = load_matrix(args.p10) if args.p10 else None
    if p10 is not None:
        save_matrix(args.out / "L_one_minus_P10.csv", score_one_minus_p10(p10))
        save_matrix(args.out / "L_hybrid_C0_P10.csv", score_hybrid(c0, p10))
        summary["p10"] = str(args.p10)
    if args.count_run2 is not None:
        u = load_matrix(args.u) if args.u else None
        f = load_matrix(args.f) if args.f else None
        ev = evaluate_leftover(
            c0,
            load_matrix(args.count_run2),
            p10=p10,
            u=u,
            f=f,
            floor=args.floor,
            cut=args.cut,
        )
        summary["evaluation"] = ev
        summary["count_run2"] = str(args.count_run2)
    _json(args.out / "instability.json", summary)
    return 0


def cmd_calibrate_leftover(args) -> int:
    """Paper 2: fit P_L=phi(L) on leftover y (never plant y). Not hat{U}."""
    args.out.mkdir(parents=True, exist_ok=True)
    if args.loso_json:
        # JSON list of {key, count_run1, count_run2}
        subjects = json.loads(Path(args.loso_json).read_text())
        loaded = []
        for s in subjects:
            loaded.append(
                {
                    "key": s["key"],
                    "c_run1": load_matrix(s["count_run1"]),
                    "c_run2": load_matrix(s["count_run2"]),
                }
            )
        rec = {"edgeuq_version": __version__} | loso_PL(loaded, floor=args.floor, cut=args.cut)
        _json(args.out / "P_L_loso.json", rec)
        return 0

    if not args.count_run1 or not args.count_run2:
        raise SystemExit("give --count-run1/--count-run2 for a single fit, or --loso-json")
    c0, c1 = load_matrix(args.count_run1), load_matrix(args.count_run2)
    L, y, _ = leftover_y(c0, c1, floor=args.floor, cut=args.cut)
    model = fit_PL(L, y)
    if args.stratum:
        model = attach_stratum_metadata(model, stratum=args.stratum, version=args.map_version)
    (args.out / "P_L_model.json").write_text(json.dumps(model, indent=2) + "\n")
    ev = evaluate_PL(L, y, model)
    # in-sample metrics are for sanity only — prefer LOSO
    rec = {
        "edgeuq_version": __version__,
        "warning": "in_sample_fit; use --loso-json for held-out P_L",
        "symbol": "P_L",
        "not_hatU": True,
        "not_plant_y": True,
        "evaluation_in_sample": ev,
        "stratum": model.get("stratum"),
    }
    _json(args.out / "P_L_fit.json", rec)
    return 0


def cmd_apply_leftover_map(args) -> int:
    """Apply a shipped site/protocol leftover map to ses-1 L. Not hat{U}."""
    args.out.mkdir(parents=True, exist_ok=True)
    model = load_leftover_map(args.map)
    c0 = load_matrix(args.count_run1)
    L = score_neg_log_count(c0)
    try:
        p = apply_leftover_map(
            L,
            model,
            stratum=args.stratum,
            force_cross_site=args.force_cross_site,
        )
    except StratumMismatchError as e:
        print(str(e), file=__import__("sys").stderr)
        return 2
    save_matrix(args.out / "P_L.csv", p)
    save_matrix(args.out / "L.csv", L)
    rec = {
        "edgeuq_version": __version__,
        "symbol": "P_L",
        "not_hatU": True,
        "map": str(args.map),
        "map_stratum": model.get("stratum"),
        "case_stratum": args.stratum,
        "force_cross_site": bool(args.force_cross_site),
        "count_run1": str(args.count_run1),
    }
    _json(args.out / "P_L_apply.json", rec)
    return 0


def cmd_demo(args) -> int:
    from .demo import write_demo

    p = write_demo(args.out / "inputs", n=args.streamlines)
    run = args.out / "run-1"
    common = ["--tck", str(p["tck"]), "--nodes", str(p["nodes"]), "--out", str(run)]
    main(["shift", *common])
    main(["jitter", *common, "--sigma", "0.5", "--samples", "8"])
    main(["instability", "--count-run1", str(run / "C0.csv"), "--out", str(run)])
    return main(["rescan", "--pair", str(run), str(p["count_a"]), str(p["count_b"]), "--out", str(args.out / "rescan.json")])


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="edgeuq", description="Edge-level uncertainty for streamline-count connectomes.")
    ap.add_argument("--version", action="version", version=f"edgeuq {__version__}")
    sub = ap.add_subparsers(
        dest="cmd",
        required=True,
        metavar="{shift,jitter,rescan,instability,calibrate-leftover,apply-leftover-map,demo}",
    )

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

    p = sub.add_parser(
        "instability",
        aliases=["leftover", "L"],
        help="leftover-instability score L=-log1p(C0) (Dir 5); separate from shift/jitter; not hat{U}",
    )
    p.add_argument("--count-run1", type=Path, required=True, help="ses-1 count matrix (CSV)")
    p.add_argument("--out", type=Path, required=True, help="output directory (writes L.csv, instability.json)")
    p.add_argument("--count-run2", type=Path, default=None, help="optional ses-2 counts to evaluate AUROC vs leftover y")
    p.add_argument("--p10", type=Path, default=None, help="optional P(count>=floor) matrix for hybrid scores")
    p.add_argument("--u", type=Path, default=None, help="optional U.csv for comparison AUROC")
    p.add_argument("--f", type=Path, default=None, help="optional F.csv for comparison AUROC")
    p.add_argument("--floor", type=float, default=COUNT_FLOOR)
    p.add_argument("--cut", type=float, default=RELCHANGE_CUT)
    p.set_defaults(func=cmd_instability)

    p = sub.add_parser(
        "calibrate-leftover",
        aliases=["PL", "p-l"],
        help="Paper 2: calibrate P_L=phi(L) on leftover y (subject LOSO); never plant y; not hat{U}",
    )
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--count-run1", type=Path, default=None, help="single-pair ses-1 counts (in-sample demo only)")
    p.add_argument("--count-run2", type=Path, default=None, help="single-pair ses-2 counts")
    p.add_argument("--loso-json", type=Path, default=None, help="JSON list of {key,count_run1,count_run2} for LOSO")
    p.add_argument("--stratum", type=str, default=None, help="stamp site/protocol id on single-fit map (e.g. kirby21_ifod2_radial4)")
    p.add_argument("--map-version", type=str, default=None, help="optional version string written into the map JSON")
    p.add_argument("--floor", type=float, default=COUNT_FLOOR)
    p.add_argument("--cut", type=float, default=RELCHANGE_CUT)
    p.set_defaults(func=cmd_calibrate_leftover)

    p = sub.add_parser(
        "apply-leftover-map",
        aliases=["apply-PL"],
        help="apply shipped site leftover map φ_L to ses-1 L; refuses stratum mismatch unless --force-cross-site",
    )
    p.add_argument("--map", type=Path, required=True, help="isotonic leftover map JSON with stratum field")
    p.add_argument("--count-run1", type=Path, required=True, help="ses-1 count matrix")
    p.add_argument("--stratum", type=str, required=True, help="declared case stratum (must match map)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--force-cross-site",
        action="store_true",
        help="override stratum check (research only; not the product default)",
    )
    p.set_defaults(func=cmd_apply_leftover_map)

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
