# Reproducing the paper results from public data

`edgeuq` starts from a tractogram and node image. The papers produced those with QSIPrep → QSIRecon (MRtrix3 single-shell SS3T, ACT-HSVS, iFOD2, 10 million streamlines) and a FreeSurfer/FastSurfer DKT parcellation (78 nodes) warped to DWI space. Any pipeline that writes a `.tck` and a node image in the same space works; results will differ in detail.

## Data

| Role | Dataset | Notes |
|---|---|---|
| Development (U calibration, C2 planted exam) | IDEAS sub-1, sub-6 | public demo sessions |
| Development | BTC glioma, OpenNeuro `ds001226`, sub-PAT20 pre-op | |
| Scan–rescan | Kirby21 (NITRC, Multi-Modal MRI Reproducibility Resource) | two sessions per subject |
| Scan–rescan | Boekel DWI test–retest | Philips data: the last volume of each run is a trace image (b > 0, zero gradient vector) and must be removed before eddy |

Clinical confirmation sessions used in the papers are not public and are not needed to run anything here.

## Steps per session

```bash
mrconvert dkt_nodes.mif nodes.nii.gz          # if your nodes are .mif (or pass .mif directly in the container)
edgeuq u  --tck streamlines.tck --nodes nodes.nii.gz --n-nodes 78 --out S
edgeuq c2 --tck streamlines.tck --nodes nodes.nii.gz --n-nodes 78 --out S --leftover-median <median_mm>
```

For a scan–rescan pair, run both on session 1 with `--sigma 0.4981`, build session 1 and session 2 count matrices with the same pipeline, then:

```bash
edgeuq retest --pair S1/ count_ses-1.csv count_ses-2.csv [--pair ...] --out retest.json
```

## Determinism

U is deterministic. C2 is bit-reproducible for a fixed seed, NumPy version and input; the papers used seed 20261019 plus a per-session index.
