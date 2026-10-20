# connectome-edge-uq

Edge-level uncertainty for streamline-count structural connectomes.

Given a tractogram (`.tck`) and a node parcellation in the same space, `edgeuq` answers, for every edge of the count connectome: *how much would this edge's count change if the streamline endpoints were placed slightly differently?* It does this with three companion methods, one subcommand each:

| Subcommand | Method | Output | Docs |
|---|---|---|---|
| `edgeuq u` | **U**: count change under planted sub-millimetre endpoint shifts, plus a frozen calibration to P(edge moves under a 1 mm shift) | `U.csv`, `Ucal.csv` | [docs/u.md](docs/u.md) |
| `edgeuq c2` | **C2**: per-endpoint soft assignment under Gaussian endpoint jitter of measured size | `F.csv` (leakage), `Esoft.csv` (soft connectome), `P_ge_floor.csv` | [docs/c2.md](docs/c2.md) |
| `edgeuq retest` | **Retest exam**: do run-1 scores predict which edges actually change between two acquisitions of the same subject? | `retest.json` | [docs/retest.md](docs/retest.md) |

The fixed protocol (cut-offs, shifts, seeds) is in [docs/protocol.md](docs/protocol.md). How to regenerate the paper numbers from public data is in [docs/reproduce.md](docs/reproduce.md).

## Install

Container (recommended):

```bash
docker pull ghcr.io/phindagijimana/connectome-edge-uq:main
apptainer build edgeuq.sif docker://ghcr.io/phindagijimana/connectome-edge-uq:main
```

Python (3.9+; needs only NumPy and NiBabel):

```bash
pip install .
```

## Quick start

```bash
# smoke test on synthetic data (about a minute)
edgeuq demo --out demo

# real data: one session
edgeuq u  --tck sub-01_streamlines.tck --nodes sub-01_nodes.nii.gz --n-nodes 78 --out sub-01
edgeuq c2 --tck sub-01_streamlines.tck --nodes sub-01_nodes.nii.gz --n-nodes 78 --out sub-01 --sigma 0.5

# retest pairs (run-1 scores vs run-1/run-2 counts); repeat --pair for more subjects
edgeuq retest --pair sub-01 sub-01_run-1_count.csv sub-01_run-2_count.csv --out retest.json
```

With Apptainer, prefix with `apptainer run -B $PWD edgeuq.sif`.

## Inputs

- **Tractogram**: MRtrix `.tck` / `.tck.gz`, coordinates in scanner mm. Only the first and last point of each streamline are read.
- **Nodes**: integer image, 0 = background, 1..N = nodes, in the same scanner space as the tractogram. NIfTI is read directly; `.mif` needs `mrconvert` (included in the container).
- **Counts** (retest only): N×N symmetric count matrices (CSV or whitespace) for run 1 and run 2, built with the same assignment rule (radial search, 4 mm), e.g. MRtrix `tck2connectome -assignment_radial_search 4`.

Assignment follows MRtrix radial search: an endpoint takes its voxel's label, otherwise the nearest labelled voxel within 4 mm, otherwise it is unassigned.

## Status

Pre-publication. The repository is private until the preprint; at that point it becomes public and gets a Zenodo DOI (see [docs/release.md](docs/release.md)).

## Licence

Apache-2.0. If you use this software, please cite it (see `CITATION.cff`) and MRtrix3, whose assignment rule it follows.
