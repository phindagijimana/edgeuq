# edgeuq

Edge-level uncertainty for streamline-count structural connectomes.

Given a tractogram (`.tck`) and a node parcellation in the same space, `edgeuq` answers, for every edge of the count connectome: *how much would this edge's count change if the streamline endpoints were placed slightly differently?* Companion methods:

| Subcommand | Method | Output | Docs |
|---|---|---|---|
| `edgeuq shift` | **Shift sensitivity (U)**: count change under planted sub-millimetre endpoint shifts, plus a frozen calibration to P(edge moves under a 1 mm shift) | `U.csv`, `Ucal.csv` | [docs/shift.md](docs/shift.md) |
| `edgeuq jitter` | **Jitter propagation**: independent Gaussian jitter of every endpoint, of measured size, propagated through assignment | `F.csv` (leakage), `Esoft.csv`, `P_ge_floor.csv` | [docs/jitter.md](docs/jitter.md) |
| `edgeuq instability` | **Leftover instability (L)**: ses-1 score \(L=-\log(1+C_0)\); separate from U/F (not \(\hat{U}\)) | `L.csv`, `instability.json` | [docs/instability.md](docs/instability.md) |
| `edgeuq rescan` | **Rescan validation**: do run-1 scores predict which edges actually change between two acquisitions of the same subject? | `rescan.json` | [docs/rescan.md](docs/rescan.md) |

What the three are for, how they relate and what they do not cover: [docs/overview.md](docs/overview.md). The fixed protocol (cut-offs, shifts, seeds) is in [docs/protocol.md](docs/protocol.md). How to regenerate the paper numbers from public data is in [docs/reproduce.md](docs/reproduce.md).

## Install

Container (recommended):

```bash
docker pull ghcr.io/phindagijimana/edgeuq:main
apptainer build edgeuq.sif docker://ghcr.io/phindagijimana/edgeuq:main
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
edgeuq shift  --tck sub-01_streamlines.tck --nodes sub-01_nodes.nii.gz --n-nodes 78 --out sub-01
edgeuq jitter --tck sub-01_streamlines.tck --nodes sub-01_nodes.nii.gz --n-nodes 78 --out sub-01 --sigma 0.5

# leftover-instability scores (Dir 5; not hat{U})
edgeuq instability --count-run1 sub-01/C0.csv --out sub-01

# rescan pairs (run-1 scores vs run-1/run-2 counts); repeat --pair for more subjects
edgeuq rescan --pair sub-01 sub-01_run-1_count.csv sub-01_run-2_count.csv --out rescan.json
```

With Apptainer, prefix with `apptainer run -B $PWD edgeuq.sif`.

## Inputs

- **Tractogram**: MRtrix `.tck` / `.tck.gz`, coordinates in scanner mm. Only the first and last point of each streamline are read.
- **Nodes**: integer image, 0 = background, 1..N = nodes, in the same scanner space as the tractogram. NIfTI is read directly; `.mif` needs `mrconvert` (included in the container).
- **Counts** (rescan only): N×N symmetric count matrices (CSV or whitespace) for run 1 and run 2, built with the same assignment rule (radial search, 4 mm), e.g. MRtrix `tck2connectome -assignment_radial_search 4`.

Assignment follows MRtrix radial search: an endpoint takes its voxel's label, otherwise the nearest labelled voxel within 4 mm, otherwise it is unassigned.

## Status

Pre-publication. The repository is private until the preprint; at that point it becomes public and gets a Zenodo DOI (see [docs/release.md](docs/release.md)).

## Licence

Apache-2.0. If you use this software, please cite it (see `CITATION.cff`) and MRtrix3, whose assignment rule it follows.
