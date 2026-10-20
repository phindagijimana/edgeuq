"""C2: per-endpoint soft assignment under Gaussian endpoint jitter of measured size."""

from __future__ import annotations

import numpy as np

from .assign import assign_endpoints, connectome, same_edge, valid_edge
from .constants import COUNT_FLOOR, MAXWELL_MEDIAN, P_HIGH, P_LOW, RADIAL_MM


def sigma_from_median(median_mm: float) -> float:
    """Isotropic 3-D Gaussian sigma whose displacement-magnitude median equals median_mm."""
    return float(median_mm) / MAXWELL_MEDIAN


def leakage(asg0: np.ndarray, stay: np.ndarray, n_samples: int, n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """F_ij = 1 - (sum of P(stay) over streamlines hard-assigned to ij) / C0_ij. Returns (F, C0)."""
    a0, b0 = asg0[:, 0], asg0[:, 1]
    ok = valid_edge(a0, b0, n_nodes)
    ia = np.minimum(a0[ok], b0[ok]) - 1
    ib = np.maximum(a0[ok], b0[ok]) - 1
    c0 = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    kept = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    np.add.at(c0, (ia, ib), 1.0)
    np.add.at(kept, (ia, ib), stay[ok].astype(np.float64) / float(n_samples))
    c0 = c0 + c0.T
    kept = kept + kept.T
    f = np.where(c0 > 0, 1.0 - kept / np.maximum(c0, 1.0), 0.0)
    return np.clip(f, 0.0, 1.0), c0


def run_c2(
    ends: np.ndarray,
    nodes: np.ndarray,
    affine: np.ndarray,
    n_nodes: int,
    sigma_mm: float,
    n_samples: int,
    seed: int,
    radius_mm: float = RADIAL_MM,
    progress=None,
) -> dict:
    asg0 = assign_endpoints(ends, nodes, affine, radius_mm)
    a0, b0 = asg0[:, 0], asg0[:, 1]
    stay = np.zeros(ends.shape[0], dtype=np.uint16)
    cs = np.zeros((n_samples, n_nodes, n_nodes), dtype=np.float32)
    rng = np.random.default_rng(seed)
    for s in range(n_samples):
        jit = ends + rng.normal(0.0, sigma_mm, size=ends.shape).astype(np.float32)
        asg = assign_endpoints(jit, nodes, affine, radius_mm)
        del jit
        stay += same_edge(a0, b0, asg[:, 0], asg[:, 1]).astype(np.uint16)
        cs[s] = connectome(asg, n_nodes)
        if progress:
            progress(s + 1, n_samples)
    f, c0 = leakage(asg0, stay, n_samples, n_nodes)
    esoft = cs.astype(np.float64).mean(axis=0)
    p_ge = (cs >= COUNT_FLOOR).mean(axis=0)
    iu = np.triu_indices(n_nodes, 1)
    ok = valid_edge(a0, b0, n_nodes)
    invented = (esoft[iu] >= COUNT_FLOOR) & (c0[iu] < COUNT_FLOOR)
    lost = (esoft[iu] < COUNT_FLOOR) & (c0[iu] >= COUNT_FLOOR)
    uncertain = (p_ge[iu] > P_LOW) & (p_ge[iu] < P_HIGH)
    return {
        "C0": c0,
        "F": f,
        "Esoft": esoft,
        "P_ge_floor": p_ge,
        "Cs": cs,
        "summary": {
            "n_streamlines": int(ends.shape[0]),
            "n_nodes": int(n_nodes),
            "sigma_mm": float(sigma_mm),
            "n_samples": int(n_samples),
            "seed": int(seed),
            "radius_mm": float(radius_mm),
            "mean_streamline_stay": float(stay[ok].mean() / n_samples) if ok.any() else None,
            "edges_hard_ge_floor": int((c0[iu] >= COUNT_FLOOR).sum()),
            "edges_soft_ge_floor": int((esoft[iu] >= COUNT_FLOOR).sum()),
            "invented_edges": int(invented.sum()),
            "lost_edges": int(lost.sum()),
            "threshold_uncertain_edges": int(uncertain.sum()),
        },
    }
