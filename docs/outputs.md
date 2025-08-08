# Outputs & Files

## Allele Trajectories CSV

Path: by `--out` (default `allele_trajectories.csv`)

Columns:
- `generation`: generation index (0 = initial)
- `snp_id`: locus index [0..L)
- `frequency`: derived allele frequency in [0,1]
- `a1`, `b1`: trait 1 additive and GxE causal effect sizes for this SNP
- `a2`, `b2`: trait 2 additive and GxE causal effect sizes for this SNP
- `s_effect`: effective selection coefficient per allele copy (average marginal effect on relative fitness, aggregated over individuals)

### Effective selection coefficient `s_effect`

We report a per-SNP, per-generation summary of the instantaneous selection pressure mediated by both traits and the environment. For locus `j` and individual `i` with environment `E_i` and phenotypes `(z1_i, z2_i)`, the log-fitness gradients under Gaussian stabilizing selection are:

`∂ log w / ∂ z1 = −(z1 − optimum1) / ω1²`, `∂ log w / ∂ z2 = −(z2 − optimum2) / ω2²`.

The per-allele marginal effect of SNP `j` on phenotype `k` is `(a_kj + b_kj · E_i)`. We average the resulting fitness gradients across individuals using the normalized reproduction weights `w_i` (sum to 1):

`s_effect[j] = Σ_i w_i · [ (∂ log w / ∂ z1)_i · (a1_j + b1_j E_i) + (∂ log w / ∂ z2)_i · (a2_j + b2_j E_i) ]`.

This yields an “effective selection coefficient per allele copy” reflecting the total selection pressure on the SNP in that generation, not just its raw phenotype effects.

## Phenotype Timeseries CSV

Path: by `--phenotype-out` (default `phenotype_timeseries.csv`)

Columns per generation:
- `generation`, `N`, `avg_fitness`
- `mean_z1`, `var_z1`, `mean_z2`, `var_z2`
- `phen_cov`, `phen_corr`
- `real_var_add1`, `real_var_gxe1`, `real_var_add2`, `real_var_gxe2`
- `rg_add_real`, `rg_gxe_real`

Notes:
- “real_var_*” values are computed from current frequencies and effect sizes.
- Phenotypic covariance/correlation are computed over individuals for that generation.
