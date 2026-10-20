# Overview: what the three methods are for

## The problem

A streamline-count connectome edge counts streamlines whose endpoints land in two nodes. Endpoints sit near region boundaries, so a small registration or segmentation error moves many of them into a neighbouring node, and the edge count changes for reasons unrelated to the brain.

## How the three relate

| Method | Question | Error model |
|---|---|---|
| [U](u.md) | Which edges are fragile? | all endpoints shifted together by a planted sub-millimetre amount |
| [C2](c2.md) | Same, under realistic error; which edges sit on the count threshold? | independent Gaussian jitter per endpoint, of measured size |
| [Retest exam](retest.md) | Do single-scan scores predict real scan–rescan change? | none (real data) |

U and C2 are two ways to score edges; the retest exam is the test of both. U and C2's leakage F tend to rank edges similarly. C2's distinct output is the threshold view: P(count ≥ 10) under jitter.

## What it is for

Error bars on connectome edges from a single scan.

- **Group studies**: down-weight or drop fragile edges before statistics, reducing noise and the multiple-comparison burden.
- **Longitudinal and clinical work**: judge whether a connection lost between two visits is plausible change or a fragile edge flipping.
- **Quality control**: see which sessions and regions depend on registration accuracy.
- **Pipeline comparison**: compare atlases and processing choices on fragility without collecting rescans.

## What it does not do

- It covers one error source: endpoint placement relative to region boundaries.
- It says nothing about whether a streamline is anatomically real; reliable is not the same as correct.
- It does not model or correct tractography randomness; in retest data that must be separated with a reseeded-tracking control.
- The soft connectome from C2 is a diagnostic, not a corrected connectome: averaging can create or remove edges at the count floor.
