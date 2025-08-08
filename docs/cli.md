# Command Line Interface

Run `wf_codex --help` for the authoritative list. Key options:

- `-N, --pop-size <N>`: number of diploid individuals
- `-L, --loci <L>`: number of biallelic loci
- `-G, --generations <G>`: generations to simulate
- `--h2-add`, `--h2-gxe`, `--h2-add2`, `--h2-gxe2`: per-trait additive and GxE heritabilities
- `--rg-add`, `--rg-gxe`: genetic correlations (additive, GxE) in [-1, 1]
- `--env-sd`: environmental SD for E ~ N(0, env_sd^2)
- `--optimum`, `--omega`: trait 1 selection optimum & width
- `--optimum2`, `--omega2`: trait 2 selection optimum & width
- `--init-freq-min`, `--init-freq-max`: initial frequency seed range (per locus)
- `-o, --out`: allele trajectory CSV path (default: allele_trajectories.csv)
- `--phenotype-out`: phenotype timeseries CSV path (default: phenotype_timeseries.csv)
- `--seed`: random seed; omitted = random
- `--threads`: worker threads (0 = auto)
- `--pop-schedule`: repeatable `gen:size` entries to change N during the run

Example:
```bash
wf_codex -N 2000 -L 1000 -G 250 --h2-add 0.2 --h2-gxe 0.1 --rg-add 0.4 \
         --optimum 0 --omega 1 --threads 8
```

