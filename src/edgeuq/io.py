"""Inputs and outputs: MRtrix .tck endpoints, node images, count matrices."""

from __future__ import annotations

import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np


def extract_endpoints(tck_path: str | Path, max_streamlines: int | None = None) -> np.ndarray:
    """First/last XYZ (scanner mm) of each streamline in a .tck or .tck.gz; shape (n, 2, 3)."""
    opener = gzip.open if str(tck_path).endswith(".gz") else open
    firsts: list[np.ndarray] = []
    lasts: list[np.ndarray] = []
    with opener(tck_path, "rb") as fh:
        data_offset: int | None = None
        while True:
            line = fh.readline()
            if not line:
                raise ValueError(f"{tck_path}: no END in tck header")
            if line.lower().startswith(b"file:"):
                parts = line.decode("ascii", "replace").split()
                if len(parts) >= 3 and parts[1] == ".":
                    data_offset = int(parts[2])
            if line.strip() == b"END":
                break
        if data_offset is not None:
            extra = data_offset - fh.tell()
            if extra > 0:
                fh.read(extra)
            elif extra < 0:
                raise ValueError(f"{tck_path}: file offset {data_offset} is before END")
        leftover = b""
        cur_first: np.ndarray | None = None
        cur_last: np.ndarray | None = None
        done = False
        while not done:
            buf = leftover + fh.read(12 * 2_000_000)
            if len(buf) < 12:
                break
            n = len(buf) // 12
            leftover = buf[n * 12 :]
            pts = np.frombuffer(buf[: n * 12], dtype="<f4").reshape(-1, 3)
            sep = np.isnan(pts[:, 0]) | np.isinf(pts[:, 0])
            start = 0
            for k in np.flatnonzero(sep).tolist() + [pts.shape[0]]:
                if k > start:
                    if cur_first is None:
                        cur_first = pts[start].copy()
                    cur_last = pts[k - 1].copy()
                if k == pts.shape[0]:
                    break
                if np.isinf(pts[k, 0]):
                    done = True
                    break
                if cur_first is not None and cur_last is not None:
                    firsts.append(cur_first)
                    lasts.append(cur_last)
                cur_first = None
                cur_last = None
                if max_streamlines is not None and len(firsts) >= max_streamlines:
                    done = True
                    break
                start = k + 1
        if cur_first is not None and cur_last is not None:
            if max_streamlines is None or len(firsts) < max_streamlines:
                firsts.append(cur_first)
                lasts.append(cur_last)
    if not firsts:
        raise ValueError(f"{tck_path}: no streamlines")
    out = np.empty((len(firsts), 2, 3), dtype=np.float32)
    out[:, 0] = np.stack(firsts)
    out[:, 1] = np.stack(lasts)
    return out


def write_tck(path: str | Path, streamlines: list[np.ndarray]) -> None:
    """Minimal float32 little-endian .tck writer (used by the synthetic demo and tests)."""
    header = [
        "mrtrix tracks",
        "datatype: Float32LE",
        f"count: {len(streamlines)}",
    ]
    body = b"".join(
        np.vstack([np.asarray(s, dtype="<f4"), np.full((1, 3), np.nan, dtype="<f4")]).tobytes()
        for s in streamlines
    ) + np.full((1, 3), np.inf, dtype="<f4").tobytes()
    text = "\n".join(header) + "\nfile: . {off}\nEND\n"
    off = len(text.format(off=0))
    while len(text.format(off=off)) != off:
        off = len(text.format(off=off))
    with open(path, "wb") as fh:
        fh.write(text.format(off=off).encode("ascii"))
        fh.write(body)


def load_nodes(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Integer node image (0 = background, 1..N = nodes) and its voxel-to-scanner affine.

    NIfTI is read directly. MRtrix .mif needs `mrconvert` on PATH.
    """
    p = Path(path)
    if p.suffix == ".mif" or p.name.endswith(".mif.gz"):
        if shutil.which("mrconvert") is None:
            raise SystemExit(f"{p}: .mif input needs mrconvert on PATH (or convert to NIfTI first)")
        with tempfile.TemporaryDirectory() as td:
            nii = Path(td) / "nodes.nii.gz"
            subprocess.run(["mrconvert", "-quiet", "-force", str(p), str(nii)], check=True)
            img = nib.load(str(nii))
            return np.rint(img.get_fdata()).astype(np.int32), np.asarray(img.affine, dtype=np.float64)
    img = nib.load(str(p))
    return np.rint(img.get_fdata()).astype(np.int32), np.asarray(img.affine, dtype=np.float64)


def load_matrix(path: str | Path) -> np.ndarray:
    text = Path(path).read_text()
    delim = "," if "," in text.splitlines()[0] else None
    return np.loadtxt(path, delimiter=delim, dtype=np.float64)


def save_matrix(path: str | Path, mat: np.ndarray, fmt: str = "%.8g") -> None:
    np.savetxt(path, np.asarray(mat, dtype=np.float64), delimiter=",", fmt=fmt)
