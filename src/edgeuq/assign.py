"""MRtrix-style radial endpoint-to-node assignment and count connectomes."""

from __future__ import annotations

import numpy as np

from .constants import RADIAL_MM


def _world_to_voxel(xyz: np.ndarray, affine: np.ndarray) -> np.ndarray:
    inv = np.linalg.inv(affine)
    hom = np.ones((xyz.shape[0], 4), dtype=np.float64)
    hom[:, :3] = xyz
    return (inv @ hom.T).T[:, :3]


def assign_endpoints(
    endpoints: np.ndarray,
    nodes: np.ndarray,
    affine: np.ndarray,
    radius_mm: float = RADIAL_MM,
    chunk: int = 250_000,
) -> np.ndarray:
    """Node label per endpoint (0 = unassigned): voxel label if non-zero, else nearest labelled voxel within radius_mm."""
    n = endpoints.shape[0]
    out = np.zeros((n, 2), dtype=np.int32)
    shape = np.asarray(nodes.shape, dtype=np.int32)
    spacing = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
    rad_vox = np.ceil(radius_mm / np.maximum(spacing, 1e-6)).astype(int)
    offsets = np.array(
        [
            (i, j, k)
            for i in range(-rad_vox[0], rad_vox[0] + 1)
            for j in range(-rad_vox[1], rad_vox[1] + 1)
            for k in range(-rad_vox[2], rad_vox[2] + 1)
        ],
        dtype=np.int32,
    )
    for end_i in (0, 1):
        for start in range(0, n, chunk):
            sl = slice(start, min(start + chunk, n))
            vox = _world_to_voxel(endpoints[sl, end_i], affine)
            nearest = np.rint(vox).astype(np.int32)
            inside = np.all((nearest >= 0) & (nearest < shape), axis=1)
            vals = np.zeros(nearest.shape[0], dtype=np.int32)
            if np.any(inside):
                ii = nearest[inside]
                vals[inside] = nodes[ii[:, 0], ii[:, 1], ii[:, 2]]
            need = vals == 0
            if not np.any(need):
                out[sl, end_i] = vals
                continue
            sel = np.where(need)[0]
            v = vox[sel]
            cand = nearest[sel][:, None, :] + offsets[None, :, :]
            valid = np.all((cand >= 0) & (cand < shape), axis=2)
            dist = np.sqrt((((cand.astype(np.float64) - v[:, None, :]) * spacing) ** 2).sum(axis=2))
            dist = np.where(valid, dist, np.inf)
            dist = np.where(dist <= radius_mm + 1e-6, dist, np.inf)
            labels = np.zeros(cand.shape[:2], dtype=np.int32)
            if np.any(valid):
                xyz = cand.reshape(-1, 3)
                ok = valid.ravel()
                picked = np.zeros(xyz.shape[0], dtype=np.int32)
                picked[ok] = nodes[xyz[ok, 0], xyz[ok, 1], xyz[ok, 2]]
                labels = picked.reshape(cand.shape[:2])
            labels = np.where(labels > 0, labels, 0)
            dist = np.where(labels > 0, dist, np.inf)
            best = np.argmin(dist, axis=1)
            rows = np.arange(labels.shape[0])
            vals[sel] = np.where(np.isfinite(dist[rows, best]), labels[rows, best], 0)
            out[sl, end_i] = vals
    return out


def valid_edge(a: np.ndarray, b: np.ndarray, n_nodes: int) -> np.ndarray:
    return (a > 0) & (b > 0) & (a != b) & (a <= n_nodes) & (b <= n_nodes)


def same_edge(a0: np.ndarray, b0: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return ((a == a0) & (b == b0)) | ((a == b0) & (b == a0))


def connectome(assignments: np.ndarray, n_nodes: int) -> np.ndarray:
    """Symmetric streamline-count matrix (no self-connections, unassigned ends dropped)."""
    mat = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    a, b = assignments[:, 0], assignments[:, 1]
    ok = valid_edge(a, b, n_nodes)
    np.add.at(mat, (np.minimum(a[ok], b[ok]) - 1, np.maximum(a[ok], b[ok]) - 1), 1.0)
    return mat + mat.T
