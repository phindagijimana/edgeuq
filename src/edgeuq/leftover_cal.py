"""Paper 2: leftover probability P_L = phi(L) with subject-held-out calibration.

Never fit on plant y. Never name results hat{U}. U/F may be ranking baselines only.
"""

from __future__ import annotations

import numpy as np

from .constants import COUNT_FLOOR, RELCHANGE_CUT
from .instability import score_neg_log_count
from .metrics import auroc, moved_mask


def leftover_y(
    c_run1: np.ndarray,
    c_run2: np.ndarray,
    *,
    floor: float = COUNT_FLOOR,
    cut: float = RELCHANGE_CUT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (L_scores, y, eligible_mask) on upper triangle."""
    c0 = np.asarray(c_run1, dtype=np.float64)
    c1 = np.asarray(c_run2, dtype=np.float64)
    iu = np.triu_indices(c0.shape[0], 1)
    elig = np.maximum(c0, c1)[iu] >= floor
    y = moved_mask(c0, c1, cut=cut, floor=floor)[iu][elig].astype(np.float64)
    L = score_neg_log_count(c0)[iu][elig]
    return L, y, elig


def isotonic_fit(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """PAVA isotonic regression; returns (x_knots, y_means) for np.interp."""
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    order = np.argsort(x)
    x, y = x[order], y[order]
    blocks = [[float(v)] for v in y]
    xs = [[float(v)] for v in x]
    i = 0
    while i < len(blocks) - 1:
        if np.mean(blocks[i]) > np.mean(blocks[i + 1]):
            blocks[i].extend(blocks[i + 1])
            xs[i].extend(xs[i + 1])
            del blocks[i + 1], xs[i + 1]
            while i > 0 and np.mean(blocks[i - 1]) > np.mean(blocks[i]):
                blocks[i - 1].extend(blocks[i])
                xs[i - 1].extend(xs[i])
                del blocks[i], xs[i]
                i -= 1
        else:
            i += 1
    xt = np.array([float(np.mean(z)) for z in xs])
    yt = np.array([float(np.mean(z)) for z in blocks])
    return xt, yt


def apply_isotonic(x: np.ndarray, xt: np.ndarray, yt: np.ndarray) -> np.ndarray:
    return np.interp(np.asarray(x, dtype=np.float64), xt, yt, left=yt[0], right=yt[-1])


def ece(p: np.ndarray, y: np.ndarray, m: int = 10) -> float:
    p = np.asarray(p, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    if p.size == 0:
        return float("nan")
    edges = np.linspace(0.0, 1.0, m + 1)
    acc, n = 0.0, float(p.size)
    for i in range(m):
        lo, hi = edges[i], edges[i + 1]
        sel = (p >= lo) & (p < hi) if i < m - 1 else (p >= lo) & (p <= hi)
        k = int(sel.sum())
        if k:
            acc += (k / n) * abs(float(y[sel].mean()) - float(p[sel].mean()))
    return float(acc)


def brier(p: np.ndarray, y: np.ndarray) -> float:
    p = np.asarray(p, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    return float(np.mean((p - y) ** 2)) if p.size else float("nan")


def fit_PL(L_train: np.ndarray, y_train: np.ndarray) -> dict:
    xt, yt = isotonic_fit(L_train, y_train)
    return {
        "method": "isotonic",
        "X_thresholds": xt.tolist(),
        "y_thresholds": yt.tolist(),
        "target": "y_leftover",
        "not_hatU": True,
        "not_plant_calibrated": True,
    }


def predict_PL(L: np.ndarray, model: dict) -> np.ndarray:
    xt = np.asarray(model["X_thresholds"], dtype=np.float64)
    yt = np.asarray(model["y_thresholds"], dtype=np.float64)
    return apply_isotonic(L, xt, yt)


def evaluate_PL(L: np.ndarray, y: np.ndarray, model: dict) -> dict:
    p = predict_PL(L, model)
    return {
        "n": int(y.size),
        "y_rate": float(y.mean()) if y.size else None,
        "mean_P_L": float(p.mean()) if p.size else None,
        "auroc_L": float(auroc(L, y)),
        "auroc_P_L": float(auroc(p, y)),
        "ECE_P_L": ece(p, y),
        "Brier_P_L": brier(p, y),
        "symbol": "P_L",
        "not_hatU": True,
    }


def loso_PL(
    subjects: list[dict],
    *,
    floor: float = COUNT_FLOOR,
    cut: float = RELCHANGE_CUT,
) -> dict:
    """Leave-one-subject-out P_L. Each dict needs c_run1, c_run2, key."""
    folds = []
    for i, te in enumerate(subjects):
        L_tr, y_tr = [], []
        for j, tr in enumerate(subjects):
            if i == j:
                continue
            L, y, _ = leftover_y(tr["c_run1"], tr["c_run2"], floor=floor, cut=cut)
            L_tr.append(L)
            y_tr.append(y)
        if not L_tr:
            continue
        model = fit_PL(np.concatenate(L_tr), np.concatenate(y_tr))
        L_te, y_te, _ = leftover_y(te["c_run1"], te["c_run2"], floor=floor, cut=cut)
        ev = evaluate_PL(L_te, y_te, model)
        folds.append({"train_keys": [s["key"] for k, s in enumerate(subjects) if k != i], "test": te["key"], **ev})
    return {
        "protocol": "LOSO_by_subject",
        "primary_score": "L=-log1p(C0)",
        "calibrated_symbol": "P_L",
        "not_hatU": True,
        "not_plant_y": True,
        "folds": folds,
        "note": "n=2 is protocol demo only; scale required for calibration claims",
    }
