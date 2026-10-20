"""U: per-edge sensitivity of the count connectome to sub-millimetre endpoint shifts."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .assign import assign_endpoints, connectome
from .constants import AXES, COUNT_FLOOR, FIT_MM, RADIAL_MM, RELCHANGE_CUT, TEST_MM
from .metrics import apply_isotonic, auroc, load_isotonic, rel_matrix


def shift_name(axis: str, t: float) -> str:
    return f"t{axis}_{'p' if t > 0 else 'm'}{abs(t):g}mm"


def shift_endpoints(ends: np.ndarray, axis: str, t: float) -> np.ndarray:
    """Move every endpoint by -t mm along one scanner axis."""
    out = np.asarray(ends, dtype=np.float32).copy()
    out[..., AXES.index(axis)] -= np.float32(t)
    return out


def shifted_counts(
    ends: np.ndarray,
    nodes: np.ndarray,
    affine: np.ndarray,
    n_nodes: int,
    radius_mm: float = RADIAL_MM,
    cache_dir: Path | None = None,
) -> dict[str, np.ndarray]:
    maps: dict[str, np.ndarray] = {}
    for axis in AXES:
        for t in (*FIT_MM, *TEST_MM):
            name = shift_name(axis, t)
            f = cache_dir / f"{name}_count.npy" if cache_dir else None
            if f is not None and f.is_file():
                maps[name] = np.load(f).astype(np.float64)
                continue
            maps[name] = connectome(assign_endpoints(shift_endpoints(ends, axis, t), nodes, affine, radius_mm), n_nodes)
            if f is not None:
                np.save(f, maps[name].astype(np.float32))
    return maps


def u_e_mask(c0: np.ndarray, fit: list[np.ndarray], test: list[np.ndarray], floor: float = COUNT_FLOOR):
    """U = max relchange over fit shifts; e = max over 1 mm test shifts; mask = any map reaches the floor."""
    u = np.zeros_like(c0)
    for c in fit:
        u = np.maximum(u, rel_matrix(c0, c, floor))
    e = np.zeros_like(c0)
    for c in test:
        e = np.maximum(e, rel_matrix(c0, c, floor))
    mx = np.maximum.reduce([c0, *fit, *test])
    return u, e, mx >= floor


def compute_u(
    ends: np.ndarray,
    nodes: np.ndarray,
    affine: np.ndarray,
    n_nodes: int,
    radius_mm: float = RADIAL_MM,
    cache_dir: Path | None = None,
) -> dict:
    c0 = connectome(assign_endpoints(ends, nodes, affine, radius_mm), n_nodes)
    maps = shifted_counts(ends, nodes, affine, n_nodes, radius_mm, cache_dir)
    fit = [maps[shift_name(a, t)] for a in AXES for t in FIT_MM]
    test = [maps[shift_name(a, t)] for a in AXES for t in TEST_MM]
    u, e, mask = u_e_mask(c0, fit, test)
    xt, yt = load_isotonic()
    ucal = apply_isotonic(u, xt, yt)
    np.fill_diagonal(ucal, 0.0)
    iu = np.triu_indices(n_nodes, 1)
    m = mask[iu]
    y = (e[iu][m] > RELCHANGE_CUT).astype(np.float64)
    return {
        "C0": c0,
        "U": u,
        "Ucal": ucal,
        "e_1mm": e,
        "mask": mask,
        "summary": {
            "n_streamlines": int(ends.shape[0]),
            "n_nodes": int(n_nodes),
            "radius_mm": float(radius_mm),
            "n_edges_scored": int(m.sum()),
            "planted_1mm_move_rate": float(y.mean()) if y.size else None,
            "planted_1mm_auroc_U": auroc(u[iu][m], y) if y.size else None,
        },
    }
