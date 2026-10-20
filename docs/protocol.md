# Frozen protocol

These values were fixed before the results they govern were computed. They live in `src/edgeuq/constants.py`; changing any of them makes outputs incomparable with the papers.

## Names

| Method | Subcommand | Earlier internal name (lock files, paper records) |
|---|---|---|
| Shift sensitivity (U, calibrated U) | `edgeuq shift` | U, experiment C |
| Jitter propagation (leakage F, soft connectome, P(count ≥ 10)) | `edgeuq jitter` | C2 |
| Rescan validation | `edgeuq rescan` | retest exam, G |

The old subcommands `u`, `c2` and `retest` still work as aliases.

## Shared

| Item | Value |
|---|---|
| Endpoint assignment | radial search, 4 mm (MRtrix `-assignment_radial_search 4` rule) |
| Edge change | `relchange(a, b) = |a - b| / max(a, b)` |
| Moved | `relchange > 0.25` and `max(a, b) >= 10` |
| Count floor | 10 streamlines |

## Shift sensitivity (U)

| Item | Value |
|---|---|
| Fit shifts | ±0.25, ±0.5 mm on each scanner axis, all endpoints, `endpoint -= t` |
| Test shifts | ±1 mm on each scanner axis |
| U | max relchange over fit shifts (zero below the floor) |
| Scored edges | any of C0 and the 18 shifted maps >= 10 |
| Calibration | isotonic, fit once on three public development sessions (IDEAS sub-1, IDEAS sub-6, BTC sub-PAT20); frozen |

## Jitter propagation

| Item | Value |
|---|---|
| Jitter | independent isotropic Gaussian per endpoint |
| sigma | measured median displacement / 1.5382; pooled scan–rescan sigma 0.4981 mm |
| Samples | S = 32 |
| Seed | 20261019 (+ session index in the papers) |
| Planted development criterion | Spearman(F, e_1mm) > 0.5 and AUROC(F, y_1mm) > 0.7 on at least 2 of 3 development and 2 of 3 confirmation sessions |
| Threshold-uncertain | 0.1 < P(count >= 10) < 0.9 |

## Rescan validation

| Item | Value |
|---|---|
| Event | moved between run 1 and run 2 (rule above) |
| Primary | AUROC of U, Ucal and F computed on run 1 only |
| Soft connectome | supported only if `moved_frac_soft < moved_frac_hard` and `invent_frac <= 0.001` on every pair |
| Calibration | Ucal may be described as calibrated for real change only if every occupied bin (n >= 20) is within ±0.15 |
| Control | reseeded run-1 tractography, to separate tracking randomness from endpoint placement |
