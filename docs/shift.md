# `edgeuq shift` — shift sensitivity (U and calibrated U)

## Question

A structural connectome edge is the number of streamlines whose two endpoints land in two nodes. Endpoints sit near grey/white boundaries, so a small registration or segmentation error can move many of them into a neighbouring node. **U** measures, per edge, how fragile the count is to such errors, using the subject's own tractogram and parcellation.

## Method

1. Assign every endpoint to a node with radial search (4 mm) and build the reference count matrix **C0**.
2. Re-assign after moving *all* endpoints by a planted shift `t` along one scanner axis (`endpoint -= t`):
   - fit shifts: ±0.25 and ±0.5 mm on x, y, z (12 maps);
   - test shifts: ±1 mm on x, y, z (6 maps).
3. For each map `C` and edge `ij`: `rel = |C0 - C| / max(C0, C)`, set to 0 unless `max(C0, C) >= 10`.
4. **U** = max of `rel` over the 12 fit maps.
   **e_1mm** = max of `rel` over the 6 test maps; the planted event is `y = e_1mm > 0.25`.
5. **mask** = edges where any of C0 and the 18 maps reaches 10 streamlines. Only these edges are scored.
6. **Ucal** = frozen isotonic map from U to P(edge moves by more than 25% under a 1 mm shift), fit once on three public development sessions and shipped in `src/edgeuq/data/isotonic_1mm.json`. It is never refit.

U is a ranking score; Ucal is a probability *for the planted 1 mm event only*. Whether it is also calibrated for real scan–rescan change is what `edgeuq rescan` tests.

## Usage

```bash
edgeuq shift --tck T.tck --nodes nodes.nii.gz --n-nodes 78 --out OUT
```

| Option | Meaning |
|---|---|
| `--n-nodes` | number of nodes; set it (e.g. 78 for DKT) so a missing top label cannot shrink the matrix |
| `--radius` | radial search radius in mm (protocol: 4) |

## Outputs (in `OUT/`)

| File | Content |
|---|---|
| `C0.csv` | reference count matrix |
| `U.csv` | U, N×N, symmetric, in [0, 1] |
| `Ucal.csv` | calibrated U, N×N |
| `e_1mm.csv` | relative change under 1 mm test shifts (self-check target) |
| `mask.csv` | 1 where the edge is scored |
| `shift.json` | provenance plus the planted 1 mm self-check AUROC of U |
| `endpoints.npy`, `shifts/` | caches, reused on re-runs |

## Cost

19 assignments of all endpoints. About 10 to 30 minutes for 10 million streamlines on one core; memory is dominated by the endpoint array (about 240 MB per 10 million streamlines).
