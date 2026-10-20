"""Synthetic inputs so the container can be smoke-tested without any real data."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from .assign import assign_endpoints, connectome
from .io import save_matrix, write_tck

SHAPE = (48, 48, 48)
VOXEL_MM = 2.0
GRID = 3


def make_nodes() -> tuple[np.ndarray, np.ndarray]:
    """GRID^3 cubic nodes separated by a 1-voxel unlabelled gap, centred in the volume."""
    nodes = np.zeros(SHAPE, dtype=np.int16)
    side = (SHAPE[0] - 8) // GRID
    label = 1
    for i in range(GRID):
        for j in range(GRID):
            for k in range(GRID):
                lo = 4 + np.array([i, j, k]) * side
                nodes[lo[0] : lo[0] + side - 1, lo[1] : lo[1] + side - 1, lo[2] : lo[2] + side - 1] = label
                label += 1
    aff = np.diag([VOXEL_MM, VOXEL_MM, VOXEL_MM, 1.0])
    aff[:3, 3] = -VOXEL_MM * (np.asarray(SHAPE) - 1) / 2.0
    return nodes, aff


def make_streamlines(nodes: np.ndarray, aff: np.ndarray, n: int, rng: np.random.Generator) -> list[np.ndarray]:
    n_nodes = int(nodes.max())
    centres = np.array([np.argwhere(nodes == k).mean(axis=0) for k in range(1, n_nodes + 1)])
    centres_mm = centres @ aff[:3, :3].T + aff[:3, 3]
    strength = rng.gamma(0.6, 1.0, size=(n_nodes, n_nodes))
    strength = np.triu(strength, 1)
    p = strength.ravel() / strength.sum()
    picks = rng.choice(p.size, size=n, p=p)
    a, b = np.unravel_index(picks, strength.shape)
    half = VOXEL_MM * ((SHAPE[0] - 8) // GRID) / 2.0
    ea = centres_mm[a] + rng.uniform(-1.15 * half, 1.15 * half, size=(n, 3))
    eb = centres_mm[b] + rng.uniform(-1.15 * half, 1.15 * half, size=(n, 3))
    mid = 0.5 * (ea + eb)
    return [np.stack([ea[i], mid[i], eb[i]]) for i in range(n)]


def write_demo(out: Path, n: int = 40_000, seed: int = 1) -> dict[str, Path]:
    """nodes.nii.gz, run-1.tck and two count matrices that differ by realistic endpoint noise."""
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    nodes, aff = make_nodes()
    nib.save(nib.Nifti1Image(nodes, aff), str(out / "nodes.nii.gz"))
    sl = make_streamlines(nodes, aff, n, rng)
    write_tck(out / "run-1.tck", sl)
    n_nodes = int(nodes.max())
    ends1 = np.stack([np.stack([s[0], s[-1]]) for s in sl]).astype(np.float32)
    c1 = connectome(assign_endpoints(ends1, nodes, aff), n_nodes)
    keep = rng.random(n) < 0.9
    ends2 = ends1[keep] + rng.normal(0.0, 0.8, size=ends1[keep].shape).astype(np.float32)
    c2 = connectome(assign_endpoints(ends2, nodes, aff), n_nodes)
    save_matrix(out / "count_run-1.csv", c1, fmt="%d")
    save_matrix(out / "count_run-2.csv", c2, fmt="%d")
    return {
        "nodes": out / "nodes.nii.gz",
        "tck": out / "run-1.tck",
        "count_a": out / "count_run-1.csv",
        "count_b": out / "count_run-2.csv",
    }
