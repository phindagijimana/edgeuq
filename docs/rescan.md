# `edgeuq rescan` — rescan validation

## Question

Shift sensitivity (U) and jitter propagation (F) are computed from a single session. Do they predict which edges *actually* change when the same subject is scanned again and the whole pipeline is re-run?

## Design

For each subject with two acquisitions:

- Run `edgeuq shift` (and optionally `edgeuq jitter`) on **run 1 only**.
- Build count matrices for run 1 and run 2 with the same pipeline and the same assignment rule.
- The event is `y_ij = relchange(C_run1, C_run2) > 0.25` among edges with `max(C_run1, C_run2) >= 10`.

The scores never see run 2.

## Reported per pair

| Field | Meaning |
|---|---|
| `auroc_U`, `auroc_Ucal`, `auroc_F` | ranking of real movers by each run-1 score |
| `auroc_lowcount` | baseline: rank by `-C_run1` (small edges move more) |
| `move_rate` | fraction of scored edges that moved |
| `moved_frac_hard` / `moved_frac_soft` | moved fraction using C_run1 vs using Esoft as the run-1 connectome |
| `invented_edges`, `invent_frac` | Esoft >= 10 where C_run1 < 10 |
| `threshold_crossing_rate` | edges that cross the floor of 10 between runs |
| `threshold_uncertain_frac`, `crossings_on_uncertain_edges` | how many edges jitter propagation flags as threshold-uncertain, and what share of real crossings they capture |

## Pooled

- `calibration_Ucal`: ECE (10 equal-width bins) of Ucal against real movers, plus per-bin mean P vs observed rate; `occupied_bins_within_band` is true if every bin with at least 20 edges is within ±0.15.
- `soft_connectome_criterion`: true if, for every pair, the soft connectome moves less than the hard one *and* invents at most 0.1% of scored edges.

Interpretation rules are fixed in advance in [protocol.md](protocol.md). Note that run-to-run change includes tractography randomness, not only endpoint placement; a reseeded run-1 tractogram is the control for that.

## Usage

```bash
edgeuq rescan \
  --pair sub-01/ sub-01_run-1_count.csv sub-01_run-2_count.csv \
  --pair sub-02/ sub-02_run-1_count.csv sub-02_run-2_count.csv \
  --out rescan.json
```

Each `RUN1_DIR` must contain `U.csv` and `Ucal.csv`; `F.csv`, `Esoft.csv` and `P_ge_floor.csv` are used if present.
