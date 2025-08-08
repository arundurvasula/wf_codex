# Outputs & Files

## Allele Trajectories CSV

Path: by `--out` (default `allele_trajectories.csv`)

Columns:
- `generation`: generation index (0 = initial)
- `snp_id`: locus index [0..L)
- `frequency`: derived allele frequency in [0,1]

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

