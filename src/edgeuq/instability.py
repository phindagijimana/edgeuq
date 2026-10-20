"""Leftover / scan-rescan instability scores (Direction 5).

Separate from shift (U) and jitter (F). Primary ses-1-only score is -log1p(C0):
weak edges are more likely to move between visits. Do not rename this as hat{U}.
"""

from __future__ import annotations

import numpy as np

from .constants import COUNT_FLOOR, RELCHANGE_CUT
from .metrics import auroc, moved_mask


def score_neg_log_count(c_run1: np.ndarray) -> np.ndarray:
    """Ses-1-only leftover-instability score: higher = weaker = more unstable."""
    c = np.asarray(c_run1, dtype=np.float64)
    return -np.log1p(np.maximum(c, 0.0))


def score_one_minus_p10(p10: np.ndarray) -> np.ndarray:
    """Optional: 1 - P(count>=10) from jitter samples (ses-1 only)."""
    return 1.0 - np.asarray(p10, dtype=np.float64)


def score_hybrid(c_run1: np.ndarray, p10: np.ndarray) -> np.ndarray:
    """Ses-1 count + threshold uncertainty (usually ~same as count alone)."""
    return score_neg_log_count(c_run1) + score_one_minus_p10(p10)


def evaluate_leftover(
    c_run1: np.ndarray,
    c_run2: np.ndarray,
    *,
    p10: np.ndarray | None = None,
    u: np.ndarray | None = None,
    f: np.ndarray | None = None,
    floor: float = COUNT_FLOOR,
    cut: float = RELCHANGE_CUT,
) -> dict:
    """AUROC of leftover-instability scores (and optional U/F) vs leftover y."""
    c0 = np.asarray(c_run1, dtype=np.float64)
    c1 = np.asarray(c_run2, dtype=np.float64)
    iu = np.triu_indices(c0.shape[0], 1)
    elig = np.maximum(c0, c1)[iu] >= floor
    y = moved_mask(c0, c1, cut=cut, floor=floor)[iu][elig].astype(np.float64)
    scores = {
        "L_neg_log1p_C0": score_neg_log_count(c0)[iu][elig],
    }
    if p10 is not None:
        p = np.asarray(p10, dtype=np.float64)
        scores["L_one_minus_P10"] = score_one_minus_p10(p)[iu][elig]
        scores["L_hybrid_C0_P10"] = score_hybrid(c0, p)[iu][elig]
    if u is not None:
        scores["U_shift"] = np.asarray(u, dtype=np.float64)[iu][elig]
    if f is not None:
        scores["F_jitter"] = np.asarray(f, dtype=np.float64)[iu][elig]

    aurocs = {k: float(auroc(v, y)) for k, v in scores.items()}
    return {
        "n_eligible": int(elig.sum()),
        "y_rate": float(y.mean()) if y.size else None,
        "count_floor": floor,
        "relchange_cut": cut,
        "auroc": aurocs,
        "primary": "L_neg_log1p_C0",
        "not_hatU": True,
        "distinct_from_shift_jitter": True,
    }
