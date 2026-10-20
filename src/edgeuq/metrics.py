"""Edge change, ranking and calibration metrics."""

from __future__ import annotations

import json
from importlib import resources

import numpy as np

from .constants import CAL_BAND, CAL_MIN_BIN_N, COUNT_FLOOR, ECE_BINS, EPS, RELCHANGE_CUT


def relchange(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    return np.abs(a - b) / np.maximum(np.maximum(a, b), EPS)


def rel_matrix(ref: np.ndarray, other: np.ndarray, floor: float = COUNT_FLOOR) -> np.ndarray:
    """relchange, zeroed where neither count reaches the floor."""
    return np.where(np.maximum(ref, other) >= floor, relchange(ref, other), 0.0)


def moved_mask(c_a: np.ndarray, c_b: np.ndarray, cut: float = RELCHANGE_CUT, floor: float = COUNT_FLOOR) -> np.ndarray:
    return (relchange(c_a, c_b) > cut) & (np.maximum(c_a, c_b) >= floor)


def auroc(score: np.ndarray, y: np.ndarray) -> float:
    """Mann-Whitney AUROC with half credit for ties."""
    s = np.asarray(score, dtype=np.float64).ravel()
    yy = np.asarray(y, dtype=np.float64).ravel()
    pos = s[yy > 0.5]
    neg = np.sort(s[yy <= 0.5])
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    lt = np.searchsorted(neg, pos, side="left")
    le = np.searchsorted(neg, pos, side="right")
    return float((lt + 0.5 * (le - lt)).mean() / neg.size)


def bin_stats(p: np.ndarray, y: np.ndarray, m: int = ECE_BINS) -> list[dict]:
    p = np.asarray(p, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    edges = np.linspace(0.0, 1.0, m + 1)
    out = []
    for i in range(m):
        lo, hi = edges[i], edges[i + 1]
        sel = (p >= lo) & ((p < hi) if i < m - 1 else (p <= hi))
        k = int(sel.sum())
        out.append(
            {
                "bin": i + 1,
                "lo": float(lo),
                "hi": float(hi),
                "n": k,
                "mean_p": None if k == 0 else float(p[sel].mean()),
                "mean_y": None if k == 0 else float(y[sel].mean()),
            }
        )
    return out


def ece(p: np.ndarray, y: np.ndarray, m: int = ECE_BINS) -> float:
    """Equal-width expected calibration error (Naeini 2015); empty bins drop."""
    p = np.asarray(p, dtype=np.float64).ravel()
    if p.size == 0:
        return float("nan")
    return float(
        sum(b["n"] / p.size * abs(b["mean_y"] - b["mean_p"]) for b in bin_stats(p, y, m) if b["n"] > 0)
    )


def occupied_bins_within_band(p: np.ndarray, y: np.ndarray) -> tuple[bool, list[dict]]:
    bins = bin_stats(p, y)
    ok = all(abs(b["mean_y"] - b["mean_p"]) <= CAL_BAND for b in bins if b["n"] >= CAL_MIN_BIN_N)
    return ok, bins


def load_isotonic() -> tuple[np.ndarray, np.ndarray]:
    """Frozen isotonic map U -> P(planted 1 mm move); fit on public development sessions only."""
    raw = json.loads(resources.files("edgeuq.data").joinpath("isotonic_1mm.json").read_text())
    return np.asarray(raw["X_thresholds"], dtype=np.float64), np.asarray(raw["y_thresholds"], dtype=np.float64)


def apply_isotonic(u: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    uu = np.clip(np.asarray(u, dtype=np.float64), float(x[0]), float(x[-1]))
    return np.interp(uu, x, y)
