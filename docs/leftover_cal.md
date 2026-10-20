# Leftover calibration (`edgeuq calibrate-leftover`) — Paper 2

Calibrate \(P_L=\phi(L)\) where \(L=-\log(1+C_0)\) and labels are **leftover** scan–rescan moves — **not** planted 1 mm events.

## Rules

- Train/test by **subject** (`--loso-json`).
- Never fit on plant \(y\).
- Never name \(P_L\) as \(\hat{U}\) or plant \(U_{\mathrm{cal}}\).

## CLI

```bash
# preferred: leave-one-subject-out
edgeuq calibrate-leftover --loso-json subjects.json --out out_PL/

# subjects.json:
# [{"key":"kirby/KKI849","count_run1":".../ses-1/dkt_connectome_count.csv",
#   "count_run2":".../ses-2/dkt_connectome_count.csv"}, ...]
```

Writes `P_L_loso.json` with per-fold AUROC/ECE/Brier.

## Status

Pilot n=2 is a **protocol demo** only. Calibration claims need scale (`scale descriptive`).
