# Leftover instability (`edgeuq instability`)

**Direction 5.** A ses-1-only score for which edges are likely to change between acquisitions. This is a **second method**, not a retune of shift (`U`) or jitter (`F`). Do **not** call it \(\hat{U}\).

## Primary score

\[
L_{ij} = -\log(1 + C^{(1)}_{ij})
\]

Weak ses-1 edges get higher \(L\) (more unstable). On the two GOLD pilots this ranked leftover movers at AUROC ≈ 0.70, beating \(U\)/\(F\) (≈ 0.58–0.60).

## Optional companions

| Score | Formula | Role |
|---|---|---|
| `L_one_minus_P10` | \(1 - P(C \ge 10)\) from jitter samples | threshold uncertainty alone (weak) |
| `L_hybrid_C0_P10` | \(L + (1-P_{10})\) | usually ≈ \(L\) alone |

## CLI

```bash
# scores only
edgeuq instability --count-run1 C0.csv --out out/

# + AUROC vs leftover y (and optional U/F comparison)
edgeuq instability --count-run1 C0.csv --count-run2 C1.csv \
  --u U.csv --f F.csv --p10 P_ge_floor.csv --out out/
```

Writes `L.csv` (alias of `L_neg_log1p_C0.csv`) and `instability.json`. Frozen floor/cut match the protocol (`count ≥ 10`, relative change ≥ 0.25).

## What this is not

- Not plant-sensitivity \(U\).
- Not a reason to retune shift/jitter thresholds.
- Not leftover ECE / GOLD calibration (those stay fail/blocked as reported).
