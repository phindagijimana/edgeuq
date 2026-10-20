"""Retest exam: do run-1 scores predict which edges change between two real acquisitions?"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .constants import COUNT_FLOOR, INVENT_FRAC, P_HIGH, P_LOW, RELCHANGE_CUT
from .io import load_matrix
from .metrics import auroc, ece, moved_mask, occupied_bins_within_band, relchange


def load_run(run_dir: Path) -> dict[str, np.ndarray | None]:
    """Read `edgeuq u` and (optionally) `edgeuq c2` outputs written to the same directory."""
    def opt(name: str) -> np.ndarray | None:
        p = run_dir / name
        return load_matrix(p) if p.is_file() else None

    out = {k: opt(f"{k}.csv") for k in ("U", "Ucal", "F", "Esoft", "P_ge_floor")}
    if out["U"] is None or out["Ucal"] is None:
        raise SystemExit(f"{run_dir}: needs U.csv and Ucal.csv from `edgeuq u`")
    return out


def pair_scores(run: dict, c_a: np.ndarray, c_b: np.ndarray, label: str) -> dict:
    n = c_a.shape[0]
    if c_b.shape != c_a.shape or run["U"].shape != c_a.shape:
        raise SystemExit(f"{label}: matrix shapes differ ({run['U'].shape}, {c_a.shape}, {c_b.shape})")
    iu = np.triu_indices(n, 1)
    m = np.maximum(c_a, c_b)[iu] >= COUNT_FLOOR
    y = moved_mask(c_a, c_b)[iu][m].astype(np.float64)
    rec: dict = {
        "pair": label,
        "n_edges": int(m.sum()),
        "move_rate": float(y.mean()) if y.size else None,
        "auroc_U": auroc(run["U"][iu][m], y),
        "auroc_Ucal": auroc(run["Ucal"][iu][m], y),
        "auroc_lowcount": auroc(-c_a[iu][m], y),
    }
    if run["F"] is not None:
        rec["auroc_F"] = auroc(run["F"][iu][m], y)
    if run["Esoft"] is not None:
        es = run["Esoft"]
        m_soft = np.maximum(es, c_b)[iu] >= COUNT_FLOOR
        invented = int(((es[iu] >= COUNT_FLOOR) & (c_a[iu] < COUNT_FLOOR)).sum())
        rec |= {
            "moved_frac_hard": rec["move_rate"],
            "moved_frac_soft": float((relchange(es, c_b)[iu][m_soft] > RELCHANGE_CUT).mean()) if m_soft.any() else None,
            "invented_edges": invented,
            "invent_frac": invented / float(max(int(m.sum()), 1)),
        }
    cross = (c_a[iu][m] >= COUNT_FLOOR) != (c_b[iu][m] >= COUNT_FLOOR)
    rec["threshold_crossing_rate"] = float(cross.mean()) if cross.size else None
    if run["P_ge_floor"] is not None:
        p = run["P_ge_floor"][iu][m]
        flag = (p > P_LOW) & (p < P_HIGH)
        rec["threshold_uncertain_frac"] = float(flag.mean()) if flag.size else None
        rec["crossings_on_uncertain_edges"] = float(flag[cross].mean()) if cross.any() else None
    rec["_cal"] = (run["Ucal"][iu][m], y)
    return rec


def exam(pairs: list[tuple[Path, Path, Path]]) -> dict:
    rows = []
    for run_dir, ca_path, cb_path in pairs:
        rows.append(pair_scores(load_run(run_dir), load_matrix(ca_path), load_matrix(cb_path), str(run_dir)))
    p = np.concatenate([r["_cal"][0] for r in rows])
    y = np.concatenate([r["_cal"][1] for r in rows])
    within, bins = occupied_bins_within_band(p, y)
    soft = [r for r in rows if r.get("moved_frac_soft") is not None]
    return {
        "pairs": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows],
        "calibration_Ucal": {
            "n_edges": int(p.size),
            "move_rate": float(y.mean()) if y.size else None,
            "ECE": ece(p, y),
            "occupied_bins_within_band": within,
            "bins": bins,
        },
        "soft_connectome_criterion": (
            all(r["moved_frac_soft"] < r["moved_frac_hard"] and r["invent_frac"] <= INVENT_FRAC for r in soft)
            if soft and len(soft) == len(rows)
            else None
        ),
    }
