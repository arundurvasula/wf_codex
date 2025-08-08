# wf_codex — Wright–Fisher Simulator with Stabilizing Selection, GxE, and Two Correlated Traits

A fast, parallel Rust simulator for diploid Wright–Fisher dynamics with:

- Two quantitative traits under Gaussian stabilizing selection
- Additive and GxE genetic architectures (environment modulates genetic effects)
- Controllable genetic correlations between traits (additive and GxE components)
- Dynamic population size schedules across generations
- CSV allele-frequency trajectories, phenotype timeseries, and live progress reporting

Built for speed, reproducibility, and clarity — ideal for experiments, teaching, and method prototyping.

---

## Features

- Parallelized compute: threads for environment sampling, phenotype computation, and reproduction
- Stabilizing selection: per-trait optima (`--optimum`, `--optimum2`) and widths (`--omega`, `--omega2`)
- Cross-trait correlations: set `--rg-add` and `--rg-gxe` in [-1, 1]
- Initialization: per‑locus frequencies from `U(init_freq_min, init_freq_max)` (defaults 0.05–0.95)
- Reproducibility: deterministic with `--seed`; per-thread seeds derived per generation
- Output: CSV trajectories with per-SNP effects and selection, and phenotype timeseries streamed to disk
- Population schedules: change N during the run with repeatable `--pop-schedule GEN:SIZE`

---

## Installation

Prerequisites:

- Rust toolchain (latest stable recommended): https://www.rust-lang.org/tools/install

Build and test:

```bash
# Clone this repo, then
cargo build
cargo test

# Optional hygiene
cargo fmt
cargo clippy -- -D warnings
```

> Note: This crate targets Rust edition `2024` (configured in `Cargo.toml`).

---

## Documentation

- In-repo docs live under `docs/`:
  - Start at `docs/index.md`
  - Usage: `docs/usage.md`
  - CLI: `docs/cli.md`
  - Model: `docs/model.md`
  - Outputs: `docs/outputs.md`
  - Build & Release: `docs/build_and_release.md`
- API docs can be built locally via `cargo doc --no-deps --open`.
- If enabled, GitHub Pages hosts the API docs built from CI.

---

## Prebuilt Binaries

You can download ready-to-run binaries from the GitHub Releases page.

- Names: `wf_codex-<platform>-<arch>` or `wf_codex-<platform>-<arch>-<cpu_tag>`
  - Examples:
    - `wf_codex-linux-x86_64` (native to CI runner CPU)
    - `wf_codex-linux-x86_64-x86-64-v3` (portable, fast on modern x86_64)
    - `wf_codex-macos-x86_64` (Intel Macs)
    - `wf_codex-macos-arm64` (Apple Silicon/M1/M2/M3)
- Choosing the right file:
  - macOS Apple Silicon: use `-macos-arm64`.
  - macOS Intel: use `-macos-x86_64`.
  - Linux x86_64: prefer `-linux-x86_64-x86-64-v3` for broad compatibility and speed; the plain `-linux-x86_64` build is tuned to the CI runner CPU and may not run on older machines.

Download and run:

```bash
# Example for macOS Apple Silicon
curl -L -o wf_codex "https://github.com/<org>/<repo>/releases/download/<tag>/wf_codex-macos-arm64"
chmod +x wf_codex
./wf_codex --help
```

Or make it available system-wide:

```bash
install -m 0755 wf_codex /usr/local/bin/wf_codex  # may require sudo
```

---

## macOS Gatekeeper (Unknown Developer)

Because these binaries are not notarized, macOS may block them with a warning like “cannot be opened because it is from an unidentified developer.” You have a few options to allow running:

- Open via Finder:
  - Right-click the binary → Open → Open. This whitelists that one binary.
- Allow from Privacy & Security:
  - First attempt to run the binary (it will be blocked).
  - Go to System Settings → Privacy & Security → “Security” section → click “Open Anyway”.
- Remove quarantine attribute (Terminal):
  - Downloads from browsers are often marked quarantined. Clear it with:
    - `xattr -d com.apple.quarantine ./wf_codex-macos-arm64`

After one of the above, run as usual:

```bash
./wf_codex --help
```

Note: If you move or re-download the file, you may need to repeat the step.

---

## Quickstart

Run a small demo with 1,000 individuals, 500 loci, and 100 generations:

```bash
cargo run -- \
  -N 1000 -L 500 -G 100 \
  --h2-add 0.20 --h2-gxe 0.10 \
  --h2-add2 0.15 --h2-gxe2 0.05 \
  --rg-add 0.5 --rg-gxe 0.2 \
  --optimum 0.0 --omega 1.0 \
  --optimum2 0.0 --omega2 1.0 \
  --threads 0 --seed 42 \
  --out allele_trajectories.csv \
  --phenotype-out phenotype_timeseries.csv
```

Change population size mid-run with a schedule (repeatable flag):

```bash
# Grow to 2000 at gen 50, then shrink to 1500 at gen 80
cargo run -- \
  -N 1000 -L 500 -G 100 \
  --pop-schedule 50:2000 \
  --pop-schedule 80:1500 \
  --threads 0 --seed 1337
```

Live progress appears on stderr, for example:

```
[00:00:12.345] gen 42/100 | avg_w=0.913201 | h2_1(add=0.200,gxe=0.100) real(0.197,0.098) | h2_2(add=0.150,gxe=0.050) real(0.148,0.052) | rg(add=0.51,gxe=0.19) | N=1000 L=500 thr=8
```

Output CSV (`allele_trajectories.csv`) includes SNP effects and effective selection coefficients:

```csv
generation,snp_id,frequency,a1,b1,a2,b2,s_effect
0,0,0.432,-0.010,0.052,0.031,-0.014,0.004
0,1,0.517,0.022,-0.006,0.012,0.003,-0.001
...
1,0,0.429,-0.010,0.052,0.031,-0.014,0.003
1,1,0.519,0.022,-0.006,0.012,0.003,-0.002
...
```

Phenotype timeseries CSV (`phenotype_timeseries.csv`) includes per-generation stats for both traits:

```csv
generation,N,avg_fitness,mean_z1,var_z1,mean_z2,var_z2,phen_cov,phen_corr,real_var_add1,real_var_gxe1,real_var_add2,real_var_gxe2,rg_add_real,rg_gxe_real
1,1000,0.913201,-0.0021,0.3567,0.0015,0.2410,-0.0003,-0.0012,0.1972,0.0983,0.1481,0.0524,0.51,0.19
...
```
Where:
- mean_z1/mean_z2: population means of phenotypes for trait 1/2.
- var_z1/var_z2: population variances of phenotypes (per generation).
- phen_cov/phen_corr: phenotypic covariance and correlation between trait 1 and 2 in the generation.
- real_var_add*/real_var_gxe*: realized additive and GxE genetic variances computed from current allele frequencies and effect sizes.
- rg_add_real/rg_gxe_real: realized genetic correlations (additive and GxE components).

### Effective selection coefficient `s_effect`

Each SNP row in the allele CSV includes `s_effect`, an average marginal effect on relative fitness per allele copy for that generation. Under Gaussian stabilizing selection on both traits, we compute log-fitness gradients with respect to each trait and weight the SNP’s marginal phenotype effects `(a + b·E)` by those gradients, averaging across individuals using the normalized reproduction weights. Formally:

`s_effect[j] = Σ_i w_i · [ −(z1_i − optimum1)/ω1² · (a1_j + b1_j E_i) + −(z2_i − optimum2)/ω2² · (a2_j + b2_j E_i) ]`.

This reflects the total selection pressure mediated through both traits and environments, not just raw effect sizes.

---

## Sample CSV

A miniature example CSV is included at `examples/sample_allele_trajectories.csv`:

```csv
generation,snp_id,frequency
0,0,0.432
0,1,0.517
0,2,0.201
1,0,0.429
1,1,0.519
1,2,0.205
```

Use it to sanity-check downstream plotting/ingestion pipelines.

---

## CLI Overview

CLI options (auto-generated from `--help`):

| Flag | Argument | Default | Description |
|------|----------|---------|-------------|
| `-N`, `--pop-size` | `<POP_SIZE>` | — | Number of diploid individuals (effective population size) |
| `-L`, `--loci` | `<LOCI>` | — | Number of biallelic loci (SNPs) |
| `-G`, `--generations` | `<GENERATIONS>` | — | Generations to simulate (not counting generation 0 output) |
| `--h2-add` | `<H2_ADD>` | `0` | Additive heritability contribution used to scale additive effect sizes |
| `--h2-gxe` | `<H2_GXE>` | `0` | GxE heritability contribution used to scale interaction effect sizes |
| `--h2-add2` | `<H2_ADD2>` | `0` | Additive heritability for second trait |
| `--h2-gxe2` | `<H2_GXE2>` | `0` | GxE heritability for second trait |
| `--rg-add` | `<RG_ADD>` | `0` | Genetic correlation between trait1 and trait2 additive effects in [-1,1] |
| `--rg-gxe` | `<RG_GXE>` | `0` | Genetic correlation between trait1 and trait2 GxE effects in [-1,1] |
| `--env-sd` | `<ENV_SD>` | `1` | Environmental standard deviation for per-individual environment E ~ N(0, env_sd^2) |
| `--optimum` | `<OPTIMUM>` | `0` | Optimum phenotype for stabilizing selection |
| `--optimum2` | `<OPTIMUM2>` | `0` | Optimum phenotype for stabilizing selection (trait 2) |
| `--omega` | `<OMEGA>` | `1` | Selection width (omega). Smaller = stronger stabilizing selection |
| `--omega2` | `<OMEGA2>` | `1` | Selection width for second trait |
| `--init-freq-min` | `<INIT_FREQ_MIN>` | `0.05` | Minimum initial allele frequency per locus (exclusive) |
| `--init-freq-max` | `<INIT_FREQ_MAX>` | `0.95` | Maximum initial allele frequency per locus (exclusive) |
| `-o`, `--out` | `<OUT>` | `allele_trajectories.csv` | Output CSV file path for allele frequency trajectories |
| `--phenotype-out` | `<FILE>` | `phenotype_timeseries.csv` | Output CSV file path for per-generation phenotype statistics |
| `--seed` | `<SEED>` | — | Random seed (optional). If not provided, a random seed is used |
| `--threads` | `<THREADS>` | `0` | Number of worker threads to use (0 = auto-detect cores) |
| `--pop-schedule` | `<POP_SCHEDULE>` | — | Population size schedule entries like `gen:size` (repeatable) |
| `--allow-fixation` | — | `false` | Allow alleles to fix or be lost (disables boundary reintroduction) |
| `-h`, `--help` | — | — | Print help |
| `-V`, `--version` | — | — | Print version |

---

## Model Sketch

- Genotypes: diploid dosages in {0,1,2}; initial p₀ per locus sampled U(a, b)
- Phenotypes: z = Σ_j (g_j − 2p₀_j) · (a_j + b_j·E), with per-individual E ~ N(0, env_sd²)
- Fitness (two traits): unnormalized Gaussian stabilizing selection around optima:
  exp(−(z₁ − opt₁)² / (2ω₁²) − (z₂ − opt₂)² / (2ω₂)²), normalized to sampling weights
- Effect sizes: per-locus (a, b) scaled to match target per-trait heritabilities; two-trait
  effects drawn with user‑specified correlations for additive and GxE components
- Drift/selection: Wright–Fisher reproduction with fitness-weighted parent sampling

Implementation notes:

- Multithreading uses the Rust standard library (`std::thread`) and channels
- Per‑thread RNG seeds derived from the global seed and generation index
- Optional boundary behavior: by default, “edge reintroduction” keeps loci away from absorbing boundaries; pass `--allow-fixation` to disable this and allow true fixation/loss

---

## Development

Useful commands:

```bash
cargo check                # fast type-check
cargo test                 # run unit/integration tests
cargo clippy -- -D warnings
cargo fmt
```

---

## Citation / Acknowledgements

If you use this simulator in a project or paper, please consider citing the repository.
Feedback and contributions are welcome!

---

## Roadmap Ideas

- Export richer summaries (e.g., phenotype stats per generation)
- Optional recombination map and multi-chromosome structure
- More selection regimes and epistasis patterns

---

Made with Rust
