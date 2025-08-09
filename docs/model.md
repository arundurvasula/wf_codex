# Model

## Overview

- Diploid Wright–Fisher reproduction with selection
- Two quantitative traits with additive and GxE components
- Per-individual environment E ~ N(0, env_sd^2)
- Per-locus effect sizes for the two traits drawn from correlated normals
- Stabilizing selection: Gaussian fitness around trait optima with widths `omega`, `omega2`

## Phenotype

For each individual with genotype dosages `g_j ∈ {0,1,2}`:

`z = Σ_j g_j · (a_j + b_j · E)`

- `a_j`: additive effect size
- `b_j`: GxE effect size
- `E`: individual environment draw

## Heritability targets and correlations

Per-locus effects are scaled to match requested heritabilities at initialization and desired cross-trait correlations:
- `h2_add`, `h2_gxe` (trait 1), `h2_add2`, `h2_gxe2` (trait 2)
- `rg_add`, `rg_gxe` are the genetic correlations for additive and GxE components

## Fitness

Two-trait stabilizing selection (unnormalized):

`exp(−(z1 − opt1)^2 / (2 ω1^2) − (z2 − opt2)^2 / (2 ω2^2))`

Weights are normalized for parent sampling.
