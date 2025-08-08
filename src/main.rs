//! Wright–Fisher simulator with stabilizing selection, GxE, and two correlated traits.
//!
//! This binary simulates a diploid population under Wright–Fisher reproduction with:
//! - Additive genetic effects and per‑individual environment E ~ N(0, env_sd^2)
//! - GxE interaction (environment modulates additive effects)
//! - Two traits whose per‑locus effect sizes are drawn from a bivariate normal with
//!   configurable genetic correlations for additive and GxE components.
//! - Stabilizing selection around trait optima with selection widths `omega` and `omega2`.
//!
//! Parallelism
//! -----------
//! - Environment sampling, phenotype computation, and reproduction are parallelized
//!   across threads using only the Rust standard library (no external threadpools).
//! - The `--threads` flag controls the parallelism; `0` auto-detects using
//!   `std::thread::available_parallelism()`.
//!
//! Usage
//! -----
//! - Build: `cargo build`
//! - Run: `cargo run -- -N 1000 -L 500 -G 100 \
//!           --h2-add 0.2 --h2-gxe 0.1 \
//!           --h2-add2 0.15 --h2-gxe2 0.05 \
//!           --rg-add 0.5 --rg-gxe 0.2 \
//!           --optimum 0.0 --omega 1.0 \
//!           --optimum2 0.0 --omega2 1.0 \
//!           --init-freq 0.5 --threads 0`
//! - Docs: `cargo doc --open`
//!
//! Output
//! ------
//! - CSV file of allele frequency trajectories with columns:
//!   - `generation`: generation index (0 is initial state)
//!   - `snp_id`: locus index in `[0, L)`
//!   - `frequency`: derived allele frequency in `[0,1]`
//!   - `a1`,`b1`: trait 1 additive and GxE causal effect sizes for this SNP
//!   - `a2`,`b2`: trait 2 additive and GxE causal effect sizes for this SNP
//!   - `s_effect`: effective selection coefficient per allele copy (average marginal effect on relative fitness, aggregated over individuals)
//! - Live progress line on stderr each generation with timestamp, average fitness,
//!   configured and realized heritabilities for both traits, and realized genetic
//!   correlations.
//!
//! Initialization
//! --------------
//! - Per-locus initial allele frequencies are sampled from a Uniform(`--init-freq-min`,
//!   `--init-freq-max`) distribution and used to generate initial dosages via
//!   Binomial(2, p0_j). Defaults are 0.05 and 0.95.
//!
//! Reproducibility
//! ---------------
//! - Pass `--seed <u64>` for deterministic runs. When omitted, a random seed is used.
//! - Threaded steps derive per‑thread seeds from the main seed and generation.
//!
use clap::Parser;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Normal, StandardNormal};
use std::collections::BTreeMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::sync::Arc;
use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

/// Command-line arguments for the `wf_codex` simulator.
///
/// Most parameters have sensible defaults, but population size (`-N`), number of loci (`-L`),
/// and generations (`-G`) are required. See examples in the crate docs.
#[derive(Parser, Debug, Clone)]
#[command(
    name = "wf_codex",
    about = "Wright–Fisher simulator with stabilizing selection and GxE",
    version
)]
struct Args {
    /// Number of diploid individuals (effective population size)
    #[arg(short = 'N', long = "pop-size")]
    pop_size: usize,

    /// Number of biallelic loci (SNPs)
    #[arg(short = 'L', long = "loci")]
    loci: usize,

    /// Number of generations to simulate (not counting generation 0 output)
    #[arg(short = 'G', long = "generations")]
    generations: usize,

    /// Additive heritability contribution used to scale additive effect sizes
    #[arg(long = "h2-add", default_value_t = 0.0)]
    h2_add: f64,

    /// GxE heritability contribution used to scale interaction effect sizes
    #[arg(long = "h2-gxe", default_value_t = 0.0)]
    h2_gxe: f64,

    /// Additive heritability for second trait
    #[arg(long = "h2-add2", default_value_t = 0.0)]
    h2_add2: f64,

    /// GxE heritability for second trait
    #[arg(long = "h2-gxe2", default_value_t = 0.0)]
    h2_gxe2: f64,

    /// Genetic correlation between trait1 and trait2 additive effects in [-1,1]
    #[arg(long = "rg-add", default_value_t = 0.0)]
    rg_add: f64,

    /// Genetic correlation between trait1 and trait2 GxE effects in [-1,1]
    #[arg(long = "rg-gxe", default_value_t = 0.0)]
    rg_gxe: f64,

    /// Environmental standard deviation for per-individual environment E ~ N(0, env_sd^2)
    #[arg(long = "env-sd", default_value_t = 1.0)]
    env_sd: f64,

    /// Optimum phenotype for stabilizing selection
    #[arg(long = "optimum", default_value_t = 0.0)]
    optimum: f64,

    /// Optimum phenotype for stabilizing selection (trait 2)
    #[arg(long = "optimum2", default_value_t = 0.0)]
    optimum2: f64,

    /// Selection width (omega). Smaller = stronger stabilizing selection.
    #[arg(long = "omega", default_value_t = 1.0)]
    omega: f64,

    /// Selection width for second trait
    #[arg(long = "omega2", default_value_t = 1.0)]
    omega2: f64,

    /// Deprecated/ignored: initial frequencies are sampled U(0.05, 0.95) per locus.
    #[arg(long = "init-freq", default_value_t = 0.5, hide = true)]
    init_freq: f64,

    /// Minimum initial allele frequency per locus (exclusive 0, default 0.05)
    #[arg(long = "init-freq-min", default_value_t = 0.05)]
    init_freq_min: f64,

    /// Maximum initial allele frequency per locus (exclusive 1, default 0.95)
    #[arg(long = "init-freq-max", default_value_t = 0.95)]
    init_freq_max: f64,

    /// Output CSV file path for allele frequency trajectories
    #[arg(short = 'o', long = "out", default_value = "allele_trajectories.csv")]
    out: String,

    /// Output CSV file path for per-generation phenotype statistics
    #[arg(long = "phenotype-out", default_value = "phenotype_timeseries.csv")]
    phenotype_out: String,

    /// Random seed (optional). If not provided, a random seed is used.
    #[arg(long = "seed")]
    seed: Option<u64>,

    /// Number of worker threads to use (0 = auto-detect cores)
    #[arg(long = "threads", default_value_t = 0)]
    threads: usize,

    /// Optional population size schedule entries like "gen:size" (e.g., 50:2000).
    /// Repeat the flag to provide multiple generation changes.
    #[arg(long = "pop-schedule")]
    pop_schedule: Vec<String>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    if !(0.0..=1.0).contains(&args.init_freq) {
        return Err("--init-freq must be in [0,1]".into());
    }
    if args.loci == 0 || args.pop_size == 0 {
        return Err("--loci and --pop-size must be > 0".into());
    }
    if args.omega <= 0.0 {
        return Err("--omega must be > 0".into());
    }
    if args.omega2 <= 0.0 {
        return Err("--omega2 must be > 0".into());
    }
    if args.env_sd < 0.0 {
        return Err("--env-sd must be >= 0".into());
    }
    if !(-1.0..=1.0).contains(&args.rg_add) {
        return Err("--rg-add must be in [-1,1]".into());
    }
    if !(-1.0..=1.0).contains(&args.rg_gxe) {
        return Err("--rg-gxe must be in [-1,1]".into());
    }
    if !(0.0..1.0).contains(&args.init_freq_min) || !(0.0..1.0).contains(&args.init_freq_max) {
        return Err("--init-freq-min and --init-freq-max must be in (0,1)".into());
    }
    if !matches!(
        args.init_freq_min.partial_cmp(&args.init_freq_max),
        Some(std::cmp::Ordering::Less)
    ) {
        return Err("--init-freq-min must be < --init-freq-max".into());
    }

    let seed = args.seed.unwrap_or_else(|| {
        // derive a seed from OS randomness when not supplied
        // Fallback to a deterministic seed if OS RNG is unavailable
        rand::thread_rng().r#gen::<u64>()
    });
    let mut rng = StdRng::seed_from_u64(seed);

    // Initialize allele frequencies and population genotypes (N x L) with dosages in {0,1,2}
    let n = args.pop_size;
    let l = args.loci;
    // Sample per-locus initial allele frequencies p0 ~ Uniform(init_freq_min, init_freq_max)
    let uni_p = Uniform::new(args.init_freq_min, args.init_freq_max);
    let p0: Vec<f64> = (0..l).map(|_| uni_p.sample(&mut rng)).collect();
    let mut genotypes: Vec<Vec<u8>> = vec![vec![0u8; l]; n];
    for row in genotypes.iter_mut() {
        for (j, val) in row.iter_mut().enumerate() {
            *val = binomial2(&mut rng, p0[j]);
        }
    }

    // Assign SNP effect sizes for two traits with desired heritabilities and genetic correlations
    let (a1, b1, a2, b2) = assign_effects_two(
        &mut rng,
        &p0,
        args.h2_add,
        args.h2_gxe,
        args.h2_add2,
        args.h2_gxe2,
        args.rg_add,
        args.rg_gxe,
        args.env_sd,
    );

    // Prepare output writer
    let file = File::create(&args.out)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "generation,snp_id,frequency,a1,b1,a2,b2,s_effect")?;

    // Prepare phenotype stats writer
    let ph_file = File::create(&args.phenotype_out)?;
    let mut ph_writer = BufWriter::new(ph_file);
    writeln!(
        ph_writer,
        "generation,N,avg_fitness,mean_z1,var_z1,mean_z2,var_z2,phen_cov,phen_corr,real_var_add1,real_var_gxe1,real_var_add2,real_var_gxe2,rg_add_real,rg_gxe_real"
    )?;

    // Output initial frequencies and selection effect (generation 0)
    let mut cur_freqs = allele_freqs(&genotypes);
    // Sample an environment for generation 0 to compute phenotypes and selection effects
    let env0 = sample_env_parallel(n, args.env_sd, seed, 0, determine_threads(args.threads))?;
    let (z1_0, z2_0) = compute_phenotypes_two_parallel(
        &genotypes,
        &a1,
        &b1,
        &a2,
        &b2,
        &p0,
        &env0,
        determine_threads(args.threads),
    );
    let fitness0 = combined_fitness(
        &z1_0,
        &z2_0,
        args.optimum,
        args.optimum2,
        args.omega,
        args.omega2,
    );
    let s_eff0 = per_locus_selection_effects(
        &z1_0,
        &z2_0,
        &env0,
        &fitness0,
        &a1,
        &b1,
        &a2,
        &b2,
        args.optimum,
        args.optimum2,
        args.omega,
        args.omega2,
    );
    for (j, &p) in cur_freqs.iter().enumerate() {
        writeln!(
            writer,
            "0,{j},{p},{},{},{},{},{}",
            a1[j], b1[j], a2[j], b2[j], s_eff0[j]
        )?;
    }

    // Determine worker threads
    let threads = determine_threads(args.threads);

    // Parse population size schedule
    let pop_sched = parse_pop_schedule(&args.pop_schedule, args.generations)?;

    // Print simulation header with parameters
    let start = Instant::now();
    eprintln!(
        "[{}] start | N={} L={} G={} h2_add={:.4} h2_gxe={:.4} env_sd={:.4} omega={:.4} init_freq_range=[{:.3},{:.3}] seed={} threads={} out={}",
        fmt_elapsed(start.elapsed()),
        n,
        l,
        args.generations,
        args.h2_add,
        args.h2_gxe,
        args.env_sd,
        args.omega,
        args.init_freq_min,
        args.init_freq_max,
        seed,
        threads,
        &args.out
    );

    // Per-generation simulation
    // env sampling handled in parallel helper
    let mut prev_line_len = 0usize;
    let mut current_n = genotypes.len();
    for generation in 1..=args.generations {
        // Simulate environment per individual (parallel)
        let env = sample_env_parallel(current_n, args.env_sd, seed, generation as u64, threads)?;

        // Compute phenotypes (both traits) and fitness
        let (phenos1, phenos2) =
            compute_phenotypes_two_parallel(&genotypes, &a1, &b1, &a2, &b2, &p0, &env, threads);
        let avg_fit = compute_avg_fitness_two(
            &phenos1,
            &phenos2,
            args.optimum,
            args.optimum2,
            args.omega,
            args.omega2,
        );
        // Per-trait sample mean and variance of phenotypes this generation
        let (mean1, var1) = mean_var(&phenos1);
        let (mean2, var2) = mean_var(&phenos2);
        let (cov12, corr12) = cov_corr(&phenos1, &phenos2, mean1, mean2, var1, var2);
        let fitness = combined_fitness(
            &phenos1,
            &phenos2,
            args.optimum,
            args.optimum2,
            args.omega,
            args.omega2,
        );
        // Per-locus effective selection coefficient for this generation
        let s_eff = per_locus_selection_effects(
            &phenos1,
            &phenos2,
            &env,
            &fitness,
            &a1,
            &b1,
            &a2,
            &b2,
            args.optimum,
            args.optimum2,
            args.omega,
            args.omega2,
        );

        // Reproduce to create next generation (N diploid offspring)
        let next_n = pop_sched.get(&generation).copied().unwrap_or(current_n);
        let genotypes_next = reproduce_parallel_to_size(
            &genotypes,
            &fitness,
            seed,
            generation as u64,
            threads,
            next_n,
        );
        let mut genotypes = genotypes_next; // shadow and move ownership
        current_n = genotypes.len();

        // Reintroduce lost or fixed alleles to maintain polymorphism at ~1/N or 1-1/N
        reintroduce_edges(&mut rng, &mut genotypes);

        // Record frequencies
        cur_freqs = allele_freqs(&genotypes);
        let (h2a1, h2g1, h2a2, h2g2, rg_add_real, rg_gxe_real) =
            realized_heritabilities_pair(&cur_freqs, &a1, &b1, &a2, &b2, args.env_sd);
        // Write phenotype stats row
        writeln!(
            ph_writer,
            "{generation},{current_n},{avg_fit},{mean1},{var1},{mean2},{var2},{cov12},{corr12},{h2a1},{h2g1},{h2a2},{h2g2},{rg_add_real},{rg_gxe_real}"
        )?;
        ph_writer.flush()?;
        for (j, &p) in cur_freqs.iter().enumerate() {
            writeln!(
                writer,
                "{generation},{j},{p},{},{},{},{},{}",
                a1[j], b1[j], a2[j], b2[j], s_eff[j]
            )?;
        }
        // Ensure data hits disk progressively for long runs
        writer.flush()?;

        // Live progress line to stderr
        let line = format!(
            "[{}] gen {}/{} | avg_w={:.6} | h2_1(add={:.3},gxe={:.3}) real({:.3},{:.3}) | h2_2(add={:.3},gxe={:.3}) real({:.3},{:.3}) | rg(add={:.2},gxe={:.2}) | N={} L={} thr={}",
            fmt_elapsed(start.elapsed()),
            generation,
            args.generations,
            avg_fit,
            args.h2_add,
            args.h2_gxe,
            h2a1,
            h2g1,
            args.h2_add2,
            args.h2_gxe2,
            h2a2,
            h2g2,
            rg_add_real,
            rg_gxe_real,
            n,
            l,
            threads
        );
        let pad = if prev_line_len > line.len() {
            prev_line_len - line.len()
        } else {
            0
        };
        eprint!("\r{}{}", line, " ".repeat(pad));
        let _ = std::io::stderr().flush();
        prev_line_len = line.len();
    }

    eprintln!(
        "Done. Wrote allele trajectories to {} (seed={}).",
        &args.out, seed
    );
    Ok(())
}

fn fmt_elapsed(d: Duration) -> String {
    let secs = d.as_secs();
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    let ms = d.subsec_millis();
    format!("{h:02}:{m:02}:{s:02}.{ms:03}")
}

/// Compute the sample mean and population variance of a slice.
fn mean_var(xs: &[f64]) -> (f64, f64) {
    let n = xs.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mut mean = 0.0f64;
    for &x in xs.iter() {
        mean += x;
    }
    mean /= n as f64;
    let mut var = 0.0f64;
    for &x in xs.iter() {
        let d = x - mean;
        var += d * d;
    }
    var /= n as f64;
    (mean, var)
}

/// Compute phenotypic covariance and correlation between two slices.
/// Uses population moments (divide by N) to match `mean_var`.
fn cov_corr(z1: &[f64], z2: &[f64], mean1: f64, mean2: f64, var1: f64, var2: f64) -> (f64, f64) {
    let n = z1.len();
    if n == 0 || n != z2.len() {
        return (0.0, 0.0);
    }
    let mut cov = 0.0f64;
    for i in 0..n {
        cov += (z1[i] - mean1) * (z2[i] - mean2);
    }
    cov /= n as f64;
    let corr = if var1 > 0.0 && var2 > 0.0 {
        cov / (var1.sqrt() * var2.sqrt())
    } else {
        0.0
    };
    (cov, corr)
}

/// Compute average combined (unnormalized) Gaussian fitness across individuals
/// for the two-trait stabilizing selection model.
fn compute_avg_fitness_two(
    phenos1: &[f64],
    phenos2: &[f64],
    optimum1: f64,
    optimum2: f64,
    omega1: f64,
    omega2: f64,
) -> f64 {
    let n = phenos1.len();
    if n == 0 {
        return 0.0;
    }
    let inv_two_om1 = 1.0 / (2.0 * omega1 * omega1);
    let inv_two_om2 = 1.0 / (2.0 * omega2 * omega2);
    let mut sum = 0.0f64;
    for i in 0..n {
        let z1 = phenos1[i];
        let z2 = phenos2[i];
        sum +=
            (-(z1 - optimum1).powi(2) * inv_two_om1 - (z2 - optimum2).powi(2) * inv_two_om2).exp();
    }
    sum / n as f64
}

/// Compute realized per-trait additive and GxE heritabilities given current frequencies.
fn realized_heritabilities_pair(
    p: &[f64],
    a1: &[f64],
    b1: &[f64],
    a2: &[f64],
    b2: &[f64],
    env_sd: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let mut var_a1 = 0.0f64;
    let mut var_a2 = 0.0f64;
    let mut cov_a = 0.0f64;
    let mut var_b1 = 0.0f64;
    let mut var_b2 = 0.0f64;
    let mut cov_b = 0.0f64;
    for j in 0..p.len() {
        let vj = 2.0 * p[j] * (1.0 - p[j]);
        var_a1 += a1[j] * a1[j] * vj;
        var_a2 += a2[j] * a2[j] * vj;
        cov_a += a1[j] * a2[j] * vj;
        var_b1 += b1[j] * b1[j] * vj;
        var_b2 += b2[j] * b2[j] * vj;
        cov_b += b1[j] * b2[j] * vj;
    }
    let h2a1 = var_a1;
    let h2a2 = var_a2;
    let h2g1 = env_sd * env_sd * var_b1;
    let h2g2 = env_sd * env_sd * var_b2;
    let rg_add = if var_a1 > 0.0 && var_a2 > 0.0 {
        cov_a / (var_a1.sqrt() * var_a2.sqrt())
    } else {
        0.0
    };
    let rg_gxe = if var_b1 > 0.0 && var_b2 > 0.0 {
        cov_b / (var_b1.sqrt() * var_b2.sqrt())
    } else {
        0.0
    };
    (h2a1, h2g1, h2a2, h2g2, rg_add, rg_gxe)
}

/// Resolve the number of worker threads to use given a requested value.
///
/// - `0` selects `std::thread::available_parallelism()`.
/// - `>=1` uses the requested value.
fn determine_threads(requested: usize) -> usize {
    if requested > 0 {
        return requested.max(1);
    }
    thread::available_parallelism()
        .unwrap_or(NonZeroUsize::new(1).unwrap())
        .get()
}

/// Parse a population size schedule from "gen:size" entries.
fn parse_pop_schedule(
    entries: &[String],
    generations: usize,
) -> Result<BTreeMap<usize, usize>, Box<dyn Error>> {
    let mut map = BTreeMap::new();
    for raw in entries {
        let (g_str, n_str) = raw
            .split_once(':')
            .ok_or_else(|| format!("invalid pop-schedule entry '{raw}', expected gen:size"))?;
        let g: usize = g_str
            .parse()
            .map_err(|_| format!("invalid generation '{g_str}' in pop-schedule"))?;
        let n: usize = n_str
            .parse()
            .map_err(|_| format!("invalid size '{n_str}' in pop-schedule"))?;
        if g == 0 || g > generations {
            return Err(
                format!("pop-schedule generation {g} out of range (1..={generations})").into(),
            );
        }
        if n == 0 {
            return Err("pop-schedule sizes must be > 0".into());
        }
        if map.insert(g, n).is_some() {
            return Err(format!("duplicate generation {g} in pop-schedule").into());
        }
    }
    Ok(map)
}

/// Draw from `Binomial(n=2, p)` returning dosages in `{0,1,2}`.
fn binomial2<R: Rng + ?Sized>(rng: &mut R, p: f64) -> u8 {
    // Returns a draw from Binomial(n=2, p)
    let p = p.clamp(0.0, 1.0);
    let mut c = 0u8;
    if rng.r#gen::<f64>() < p {
        c += 1;
    }
    if rng.r#gen::<f64>() < p {
        c += 1;
    }
    c
}

/// Draw per-locus effect sizes for a single trait to match target heritabilities.
///
/// The raw per-locus effects are sampled from `N(0,1)` and then rescaled so that
/// the additive and GxE components contribute the requested heritabilities at the
/// initial allele frequencies `p0`.
#[cfg(test)]
fn assign_effects<R: Rng + ?Sized>(
    rng: &mut R,
    p0: &[f64],
    h2_add: f64,
    h2_gxe: f64,
    env_sd: f64,
) -> (Vec<f64>, Vec<f64>) {
    let l = p0.len();
    let normal = Normal::new(0.0, 1.0).unwrap();
    let mut a0 = vec![0.0; l];
    let mut b0 = vec![0.0; l];
    for j in 0..l {
        a0[j] = normal.sample(rng);
        b0[j] = normal.sample(rng);
    }

    let mut s_add = 0.0f64;
    let mut s_gxe = 0.0f64;
    for j in 0..l {
        let vj = 2.0 * p0[j] * (1.0 - p0[j]); // Var(Binomial(2,p))
        s_add += a0[j] * a0[j] * vj;
        s_gxe += b0[j] * b0[j] * vj;
    }
    let scale_add = if h2_add > 0.0 && s_add > 0.0 {
        (h2_add / s_add).sqrt()
    } else {
        0.0
    };
    // Var(E * sum b g) = Var(E) * Var(sum b g) with zero-mean E
    let scale_gxe = if h2_gxe > 0.0 && s_gxe > 0.0 && env_sd > 0.0 {
        (h2_gxe / (s_gxe * env_sd * env_sd)).sqrt()
    } else {
        0.0
    };

    let a: Vec<f64> = a0.into_iter().map(|x| x * scale_add).collect();
    let b: Vec<f64> = b0.into_iter().map(|x| x * scale_gxe).collect();
    (a, b)
}

#[allow(clippy::too_many_arguments)]
fn assign_effects_two<R: Rng + ?Sized>(
    rng: &mut R,
    p0: &[f64],
    h2_add1: f64,
    h2_gxe1: f64,
    h2_add2: f64,
    h2_gxe2: f64,
    rg_add: f64,
    rg_gxe: f64,
    env_sd: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let l = p0.len();
    let rho_a = rg_add.clamp(-1.0, 1.0);
    let rho_b = rg_gxe.clamp(-1.0, 1.0);

    // Raw correlated normals for additive effects
    let mut a1_raw = vec![0.0; l];
    let mut a2_raw = vec![0.0; l];
    for j in 0..l {
        let u1: f64 = StandardNormal.sample(rng);
        let u2: f64 = StandardNormal.sample(rng);
        a1_raw[j] = u1;
        a2_raw[j] = rho_a * u1 + (1.0 - rho_a * rho_a).max(0.0).sqrt() * u2;
    }

    // Raw correlated normals for GxE effects
    let mut b1_raw = vec![0.0; l];
    let mut b2_raw = vec![0.0; l];
    for j in 0..l {
        let u1: f64 = StandardNormal.sample(rng);
        let u2: f64 = StandardNormal.sample(rng);
        b1_raw[j] = u1;
        b2_raw[j] = rho_b * u1 + (1.0 - rho_b * rho_b).max(0.0).sqrt() * u2;
    }

    // Scale to target heritabilities using initial frequencies
    let mut s_a1 = 0.0f64;
    let mut s_a2 = 0.0f64;
    let mut s_b1 = 0.0f64;
    let mut s_b2 = 0.0f64;
    for j in 0..l {
        let vj = 2.0 * p0[j] * (1.0 - p0[j]);
        s_a1 += a1_raw[j] * a1_raw[j] * vj;
        s_a2 += a2_raw[j] * a2_raw[j] * vj;
        s_b1 += b1_raw[j] * b1_raw[j] * vj;
        s_b2 += b2_raw[j] * b2_raw[j] * vj;
    }
    let c_a1 = if h2_add1 > 0.0 && s_a1 > 0.0 {
        (h2_add1 / s_a1).sqrt()
    } else {
        0.0
    };
    let c_a2 = if h2_add2 > 0.0 && s_a2 > 0.0 {
        (h2_add2 / s_a2).sqrt()
    } else {
        0.0
    };
    let c_b1 = if h2_gxe1 > 0.0 && s_b1 > 0.0 && env_sd > 0.0 {
        (h2_gxe1 / (s_b1 * env_sd * env_sd)).sqrt()
    } else {
        0.0
    };
    let c_b2 = if h2_gxe2 > 0.0 && s_b2 > 0.0 && env_sd > 0.0 {
        (h2_gxe2 / (s_b2 * env_sd * env_sd)).sqrt()
    } else {
        0.0
    };

    let a1: Vec<f64> = a1_raw.into_iter().map(|x| x * c_a1).collect();
    let a2: Vec<f64> = a2_raw.into_iter().map(|x| x * c_a2).collect();
    let b1: Vec<f64> = b1_raw.into_iter().map(|x| x * c_b1).collect();
    let b2: Vec<f64> = b2_raw.into_iter().map(|x| x * c_b2).collect();
    (a1, b1, a2, b2)
}

/// Compute a single phenotype for all individuals.
/// Compute a single phenotype for all individuals.
fn compute_phenotypes(
    genotypes: &[Vec<u8>],
    a: &[f64],
    b: &[f64],
    p0: &[f64],
    env: &[f64],
) -> Vec<f64> {
    let n = genotypes.len();
    let l = a.len();
    let mut pheno = vec![0.0f64; n];
    for (i, g_row) in genotypes.iter().enumerate() {
        let e = env[i];
        let mut z = 0.0f64;
        for j in 0..l {
            let g = g_row[j] as f64;
            let dev = g - 2.0 * p0[j]; // center at initial freq
            z += dev * (a[j] + b[j] * e);
        }
        pheno[i] = z;
    }
    pheno
}

//

#[allow(clippy::too_many_arguments)]
fn compute_phenotypes_two_parallel(
    genotypes: &[Vec<u8>],
    a1: &[f64],
    b1: &[f64],
    a2: &[f64],
    b2: &[f64],
    p0: &[f64],
    env: &[f64],
    threads: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n = genotypes.len();
    if n == 0 {
        return (vec![], vec![]);
    }
    if threads <= 1 {
        let z1 = compute_phenotypes(genotypes, a1, b1, p0, env);
        let z2 = compute_phenotypes(genotypes, a2, b2, p0, env);
        return (z1, z2);
    }
    let a1 = Arc::new(a1.to_vec());
    let b1 = Arc::new(b1.to_vec());
    let a2 = Arc::new(a2.to_vec());
    let b2 = Arc::new(b2.to_vec());
    let p0 = Arc::new(p0.to_vec());
    let env = Arc::new(env.to_vec());
    let genotypes = Arc::new(genotypes.to_owned());

    let (tx, rx) = mpsc::channel();
    let chunk = n.div_ceil(threads);
    for t in 0..threads {
        let start = t * chunk;
        if start >= n {
            break;
        }
        let end = ((t + 1) * chunk).min(n);
        let tx = tx.clone();
        let a1 = Arc::clone(&a1);
        let b1 = Arc::clone(&b1);
        let a2 = Arc::clone(&a2);
        let b2 = Arc::clone(&b2);
        let p0 = Arc::clone(&p0);
        let env = Arc::clone(&env);
        let genos = Arc::clone(&genotypes);
        thread::spawn(move || {
            let mut out1 = vec![0.0f64; end - start];
            let mut out2 = vec![0.0f64; end - start];
            let l = a1.len();
            for (local_i, i) in (start..end).enumerate() {
                let e = env[i];
                let g_row = &genos[i];
                let mut z1 = 0.0f64;
                let mut z2 = 0.0f64;
                for j in 0..l {
                    let g = g_row[j] as f64;
                    let dev = g - 2.0 * p0[j];
                    z1 += dev * (a1[j] + b1[j] * e);
                    z2 += dev * (a2[j] + b2[j] * e);
                }
                out1[local_i] = z1;
                out2[local_i] = z2;
            }
            let _ = tx.send((start, out1, out2));
        });
    }
    drop(tx);
    let mut z1 = vec![0.0f64; n];
    let mut z2 = vec![0.0f64; n];
    for (start, part1, part2) in rx.iter() {
        z1[start..start + part1.len()].copy_from_slice(&part1);
        z2[start..start + part2.len()].copy_from_slice(&part2);
    }
    (z1, z2)
}

/// Compute normalized Gaussian stabilizing selection weights for a single trait.
#[cfg(test)]
fn gaussian_fitness(phenotypes: &[f64], optimum: f64, omega: f64) -> Vec<f64> {
    let two_omega2 = 2.0 * omega * omega;
    // Use log-space stabilization to avoid underflow
    let mut logs: Vec<f64> = Vec::with_capacity(phenotypes.len());
    let mut max_log = f64::NEG_INFINITY;
    for &z in phenotypes.iter() {
        let s = -((z - optimum) * (z - optimum)) / two_omega2;
        logs.push(s);
        if s > max_log {
            max_log = s;
        }
    }
    let mut weights = Vec::with_capacity(phenotypes.len());
    let mut sum = 0.0f64;
    for &s in logs.iter() {
        let w = (s - max_log).exp();
        weights.push(w);
        sum += w;
    }
    if sum <= 0.0 {
        // Fallback to equal weights if numerical issues arise
        return vec![1.0; phenotypes.len()];
    }
    for w in &mut weights {
        *w /= sum;
    }
    weights
}

/// Compute normalized combined Gaussian fitness for two traits under stabilizing selection.
fn combined_fitness(
    phenos1: &[f64],
    phenos2: &[f64],
    optimum1: f64,
    optimum2: f64,
    omega1: f64,
    omega2: f64,
) -> Vec<f64> {
    let n = phenos1.len();
    assert_eq!(n, phenos2.len());
    let inv_two_om1 = 1.0 / (2.0 * omega1 * omega1);
    let inv_two_om2 = 1.0 / (2.0 * omega2 * omega2);
    let mut logs: Vec<f64> = Vec::with_capacity(n);
    let mut max_log = f64::NEG_INFINITY;
    for i in 0..n {
        let z1 = phenos1[i];
        let z2 = phenos2[i];
        let s = -((z1 - optimum1) * (z1 - optimum1)) * inv_two_om1
            - ((z2 - optimum2) * (z2 - optimum2)) * inv_two_om2;
        logs.push(s);
        if s > max_log {
            max_log = s;
        }
    }
    let mut w = vec![0.0f64; n];
    let mut sum = 0.0f64;
    for (i, wi) in w.iter_mut().enumerate() {
        let val = (logs[i] - max_log).exp();
        *wi = val;
        sum += val;
    }
    if sum <= 0.0 {
        return vec![1.0 / n as f64; n];
    }
    for wi in &mut w {
        *wi /= sum;
    }
    w
}

/// Compute per-locus effective selection coefficients (average marginal effect on relative fitness).
///
/// For each SNP j, we compute:
///   s_j = sum_i w_i_norm * [ dlogw/dz1_i * (a1_j + b1_j * E_i) + dlogw/dz2_i * (a2_j + b2_j * E_i) ]
/// where w_i_norm are the normalized reproduction weights (sum to 1 across individuals),
/// dlogw/dz_k = -(z_k - optimum_k) / omega_k^2, and E_i is the environment.
#[allow(clippy::too_many_arguments)]
fn per_locus_selection_effects(
    phenos1: &[f64],
    phenos2: &[f64],
    env: &[f64],
    fitness_norm: &[f64],
    a1: &[f64],
    b1: &[f64],
    a2: &[f64],
    b2: &[f64],
    optimum1: f64,
    optimum2: f64,
    omega1: f64,
    omega2: f64,
) -> Vec<f64> {
    let n = phenos1.len();
    let l = a1.len();
    assert_eq!(phenos2.len(), n);
    assert_eq!(env.len(), n);
    assert_eq!(fitness_norm.len(), n);
    assert_eq!(b1.len(), l);
    assert_eq!(a2.len(), l);
    assert_eq!(b2.len(), l);

    let inv_om1_sq = 1.0 / (omega1 * omega1);
    let inv_om2_sq = 1.0 / (omega2 * omega2);
    let mut dlogw1 = vec![0.0f64; n];
    let mut dlogw2 = vec![0.0f64; n];
    for i in 0..n {
        dlogw1[i] = -((phenos1[i] - optimum1) * inv_om1_sq);
        dlogw2[i] = -((phenos2[i] - optimum2) * inv_om2_sq);
    }

    let mut out = vec![0.0f64; l];
    for j in 0..l {
        let a1j = a1[j];
        let b1j = b1[j];
        let a2j = a2[j];
        let b2j = b2[j];
        let mut s = 0.0f64;
        for i in 0..n {
            let dz1_dg = a1j + b1j * env[i];
            let dz2_dg = a2j + b2j * env[i];
            let dlogw_dg = dlogw1[i] * dz1_dg + dlogw2[i] * dz2_dg;
            s += fitness_norm[i] * dlogw_dg;
        }
        out[j] = s;
    }
    out
}

/// Produce the next generation of genotypes via Wright–Fisher sampling.
#[cfg(test)]
fn reproduce<R: Rng + ?Sized>(
    rng: &mut R,
    genotypes: &[Vec<u8>],
    fitness_weights: &[f64],
) -> Vec<Vec<u8>> {
    let n = genotypes.len();
    let l = genotypes[0].len();
    // Build cumulative distribution function for parent sampling
    let mut cdf = vec![0.0f64; n];
    let mut acc = 0.0;
    for (i, c) in cdf.iter_mut().enumerate() {
        acc += fitness_weights[i];
        *c = acc;
    }
    if acc <= 0.0 {
        for (i, c) in cdf.iter_mut().enumerate() {
            *c = (i as f64 + 1.0) / n as f64;
        }
    } else {
        for c in &mut cdf {
            *c /= acc;
        }
    }

    let mut next = vec![vec![0u8; l]; n];
    let uniform01 = Uniform::new(0.0f64, 1.0f64);

    for child in next.iter_mut() {
        let p1 = sample_index(rng, &cdf, &uniform01);
        let p2 = sample_index(rng, &cdf, &uniform01);
        let g1 = &genotypes[p1];
        let g2 = &genotypes[p2];
        for j in 0..l {
            let allele1 = draw_gamete_allele(rng, g1[j]);
            let allele2 = draw_gamete_allele(rng, g2[j]);
            child[j] = (allele1 + allele2) as u8;
        }
    }

    next
}

/// Parallel version of `reproduce` using per-thread RNGs.
#[cfg(test)]
fn reproduce_parallel(
    genotypes: &[Vec<u8>],
    fitness_weights: &[f64],
    seed: u64,
    generation: u64,
    threads: usize,
) -> Vec<Vec<u8>> {
    if threads <= 1 {
        let mut rng = StdRng::seed_from_u64(seed ^ generation.wrapping_mul(0x9E37_79B9));
        return reproduce(&mut rng, genotypes, fitness_weights);
    }

    let n = genotypes.len();
    if n == 0 {
        return vec![];
    }
    let l = genotypes[0].len();

    // Build normalized CDF once and share
    let mut cdf = vec![0.0f64; n];
    let mut acc = 0.0;
    for (i, c) in cdf.iter_mut().enumerate() {
        acc += fitness_weights[i];
        *c = acc;
    }
    if acc <= 0.0 {
        for (i, c) in cdf.iter_mut().enumerate() {
            *c = (i as f64 + 1.0) / n as f64;
        }
    } else {
        for c in &mut cdf {
            *c /= acc;
        }
    }

    let cdf = Arc::new(cdf);
    let genotypes = Arc::new(genotypes.to_owned());

    let (tx, rx) = mpsc::channel();
    let chunk = n.div_ceil(threads);
    for t in 0..threads {
        let start = t * chunk;
        if start >= n {
            break;
        }
        let end = ((t + 1) * chunk).min(n);
        let tx = tx.clone();
        let cdf = Arc::clone(&cdf);
        let genos = Arc::clone(&genotypes);
        let local_seed =
            seed ^ generation.wrapping_mul(0x9E37_79B9) ^ (t as u64).wrapping_mul(0x85EB_CA6B);
        thread::spawn(move || {
            let mut rng = StdRng::seed_from_u64(local_seed);
            let uniform01 = Uniform::new(0.0f64, 1.0f64);
            let mut part = vec![vec![0u8; l]; end - start];
            for (local_i, _) in (start..end).enumerate() {
                let p1 = sample_index(&mut rng, &cdf, &uniform01);
                let p2 = sample_index(&mut rng, &cdf, &uniform01);
                let g1 = &genos[p1];
                let g2 = &genos[p2];
                let child = &mut part[local_i];
                for j in 0..l {
                    let allele1 = draw_gamete_allele(&mut rng, g1[j]);
                    let allele2 = draw_gamete_allele(&mut rng, g2[j]);
                    child[j] = (allele1 + allele2) as u8;
                }
            }
            let _ = tx.send((start, part));
        });
    }
    drop(tx);

    let mut next = vec![vec![0u8; l]; n];
    for (start, part) in rx.iter() {
        for (offset, row) in part.into_iter().enumerate() {
            next[start + offset] = row;
        }
    }
    next
}

/// Parallel reproduction producing exactly `next_n` offspring regardless of current N.
fn reproduce_parallel_to_size(
    genotypes: &[Vec<u8>],
    fitness_weights: &[f64],
    seed: u64,
    generation: u64,
    threads: usize,
    next_n: usize,
) -> Vec<Vec<u8>> {
    let parent_n = genotypes.len();
    if parent_n == 0 || next_n == 0 {
        return vec![];
    }
    let l = genotypes[0].len();

    if threads <= 1 {
        // Build CDF
        let mut cdf = vec![0.0f64; parent_n];
        let mut acc = 0.0;
        for (i, c) in cdf.iter_mut().enumerate() {
            acc += fitness_weights[i];
            *c = acc;
        }
        if acc <= 0.0 {
            for (i, c) in cdf.iter_mut().enumerate() {
                *c = (i as f64 + 1.0) / parent_n as f64;
            }
        } else {
            for c in &mut cdf {
                *c /= acc;
            }
        }
        let mut rng = StdRng::seed_from_u64(seed ^ generation.wrapping_mul(0x9E37_79B9));
        let uniform01 = Uniform::new(0.0f64, 1.0f64);
        let mut next = vec![vec![0u8; l]; next_n];
        for child in next.iter_mut() {
            let p1 = sample_index(&mut rng, &cdf, &uniform01);
            let p2 = sample_index(&mut rng, &cdf, &uniform01);
            let g1 = &genotypes[p1];
            let g2 = &genotypes[p2];
            for j in 0..l {
                let a1 = draw_gamete_allele(&mut rng, g1[j]);
                let a2 = draw_gamete_allele(&mut rng, g2[j]);
                child[j] = (a1 + a2) as u8;
            }
        }
        return next;
    }

    // Build normalized CDF once and share among threads
    let mut cdf = vec![0.0f64; parent_n];
    let mut acc = 0.0;
    for (i, c) in cdf.iter_mut().enumerate() {
        acc += fitness_weights[i];
        *c = acc;
    }
    if acc <= 0.0 {
        for (i, c) in cdf.iter_mut().enumerate() {
            *c = (i as f64 + 1.0) / parent_n as f64;
        }
    } else {
        for c in &mut cdf {
            *c /= acc;
        }
    }

    let cdf = Arc::new(cdf);
    let genotypes = Arc::new(genotypes.to_owned());

    let (tx, rx) = mpsc::channel();
    let chunk = next_n.div_ceil(threads);
    for t in 0..threads {
        let start = t * chunk;
        if start >= next_n {
            break;
        }
        let end = ((t + 1) * chunk).min(next_n);
        let tx = tx.clone();
        let cdf = Arc::clone(&cdf);
        let genos = Arc::clone(&genotypes);
        let local_seed =
            seed ^ generation.wrapping_mul(0x9E37_79B9) ^ (t as u64).wrapping_mul(0x85EB_CA6B);
        thread::spawn(move || {
            let mut rng = StdRng::seed_from_u64(local_seed);
            let uniform01 = Uniform::new(0.0f64, 1.0f64);
            let mut part = vec![vec![0u8; l]; end - start];
            for (local_i, _) in (start..end).enumerate() {
                let p1 = sample_index(&mut rng, &cdf, &uniform01);
                let p2 = sample_index(&mut rng, &cdf, &uniform01);
                let g1 = &genos[p1];
                let g2 = &genos[p2];
                let child = &mut part[local_i];
                for j in 0..l {
                    let a1 = draw_gamete_allele(&mut rng, g1[j]);
                    let a2 = draw_gamete_allele(&mut rng, g2[j]);
                    child[j] = (a1 + a2) as u8;
                }
            }
            let _ = tx.send((start, part));
        });
    }
    drop(tx);

    let mut next = vec![vec![0u8; l]; next_n];
    for (start, part) in rx.iter() {
        for (offset, row) in part.into_iter().enumerate() {
            next[start + offset] = row;
        }
    }
    next
}

/// Sample per-individual environments `E ~ N(0, env_sd^2)` in parallel.
fn sample_env_parallel(
    n: usize,
    env_sd: f64,
    seed: u64,
    generation: u64,
    threads: usize,
) -> Result<Vec<f64>, Box<dyn Error>> {
    if n == 0 {
        return Ok(vec![]);
    }
    if threads <= 1 {
        let mut rng = StdRng::seed_from_u64(seed ^ generation.wrapping_mul(0xD1B5_4A32));
        let normal_env = Normal::new(0.0, env_sd.max(0.0))?;
        let mut env = Vec::with_capacity(n);
        for _ in 0..n {
            env.push(normal_env.sample(&mut rng));
        }
        return Ok(env);
    }

    let (tx, rx) = mpsc::channel();
    let chunk = n.div_ceil(threads);
    for t in 0..threads {
        let start = t * chunk;
        if start >= n {
            break;
        }
        let end = ((t + 1) * chunk).min(n);
        let tx = tx.clone();
        let local_seed =
            seed ^ generation.wrapping_mul(0xD1B5_4A32) ^ (t as u64).wrapping_mul(0xC2B2_AE35);
        let sd = env_sd;
        thread::spawn(move || {
            let mut rng = StdRng::seed_from_u64(local_seed);
            let normal_env = Normal::new(0.0, sd.max(0.0)).unwrap();
            let mut out = Vec::with_capacity(end - start);
            for _ in start..end {
                out.push(normal_env.sample(&mut rng));
            }
            let _ = tx.send((start, out));
        });
    }
    drop(tx);

    let mut env = vec![0.0f64; n];
    for (start, part) in rx.iter() {
        env[start..start + part.len()].copy_from_slice(&part);
    }
    Ok(env)
}

/// Sample an index from a normalized CDF using inverse transform sampling.
fn sample_index<R: Rng + ?Sized>(rng: &mut R, cdf: &[f64], uniform01: &Uniform<f64>) -> usize {
    let u = uniform01.sample(rng);
    // binary search
    match cdf.binary_search_by(|x| x.partial_cmp(&u).unwrap()) {
        Ok(idx) => idx,
        Err(idx) => idx,
    }
}

/// Draw a single gamete allele from a diploid genotype dosage in `{0,1,2}`.
fn draw_gamete_allele<R: Rng + ?Sized>(rng: &mut R, genotype: u8) -> u8 {
    match genotype {
        0 => 0,
        1 => {
            if rng.r#gen::<f64>() < 0.5 {
                1
            } else {
                0
            }
        }
        2 => 1,
        _ => unreachable!(),
    }
}

/// Compute derived allele frequencies across loci for the current population genotypes.
fn allele_freqs(genotypes: &[Vec<u8>]) -> Vec<f64> {
    let n = genotypes.len();
    let l = genotypes[0].len();
    let mut counts = vec![0u32; l];
    for row in genotypes.iter() {
        for (j, &g) in row.iter().enumerate() {
            counts[j] += g as u32;
        }
    }
    counts
        .into_iter()
        .map(|c| c as f64 / (2.0 * n as f64))
        .collect()
}

/// Maintain polymorphism by reintroducing/extirpating two allele copies at fixation/loss edges.
///
/// This keeps minor allele frequency around `1/N` away from absorbing boundaries, which is
/// useful for long simulations when tracking trajectories rather than true absorbing states.
fn reintroduce_edges<R: Rng + ?Sized>(rng: &mut R, genotypes: &mut [Vec<u8>]) {
    let n = genotypes.len();
    let l = genotypes[0].len();

    // Count totals per locus
    let mut counts = vec![0i32; l];
    for row in genotypes.iter() {
        for (j, &g) in row.iter().enumerate() {
            counts[j] += g as i32;
        }
    }

    // For each locus, if count is 0 -> set to 2 (freq = 1/N). If 2N -> set to 2N-2 (freq = 1 - 1/N).
    let rng_idx = Uniform::new(0usize, n);
    for (j, total) in counts.iter().enumerate() {
        let target_total;
        if *total == 0 {
            target_total = 2i32; // introduce two copies
        } else if *total == (2 * n) as i32 {
            target_total = (2 * n) as i32 - 2; // remove two copies
        } else {
            continue;
        }

        let mut current = *total;
        // Adjust by randomly adding/removing one allele at a time
        while current != target_total {
            let i = rng_idx.sample(rng);
            let g = genotypes[i][j] as i32;
            if target_total > current {
                // add allele if possible
                if g < 2 {
                    genotypes[i][j] = (g + 1) as u8;
                    current += 1;
                }
            } else {
                // remove allele if possible
                if g > 0 {
                    genotypes[i][j] = (g - 1) as u8;
                    current -= 1;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "heavy-tests")]
    use rand::Rng as _;

    // --- Neutral Wright–Fisher helpers for tests ---
    #[cfg(feature = "heavy-tests")]
    fn neutral_fix_or_loss_time<R: Rng + ?Sized>(rng: &mut R, n: usize, p0: f64) -> (u64, bool) {
        // Simulate a single neutral biallelic locus under classical WF model
        // until absorption (fixation or loss). Returns (generations, fixed?).
        let two_n = 2 * n as u32;
        let mut count: u32 = ((p0 * two_n as f64).round() as i64).clamp(0, two_n as i64) as u32;
        let mut t: u64 = 0;
        // Simple Bernoulli-sum binomial sampler; n is modest in tests.
        let mut binomial = |nn: u32, p: f64| -> u32 {
            let p = p.clamp(0.0, 1.0);
            let mut c = 0u32;
            for _ in 0..nn {
                if rng.r#gen::<f64>() < p {
                    c += 1;
                }
            }
            c
        };
        while count > 0 && count < two_n {
            let p = count as f64 / two_n as f64;
            count = binomial(two_n, p);
            t += 1;
        }
        (t, count == two_n)
    }

    #[test]
    fn test_binomial2_bounds() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..100 {
            assert_eq!(binomial2(&mut rng, 0.0), 0);
            assert_eq!(binomial2(&mut rng, 1.0), 2);
        }
    }

    #[test]
    fn test_assign_effects_scaling() {
        let mut rng = StdRng::seed_from_u64(2);
        let l = 100usize;
        let p0 = vec![0.3f64; l];
        let h2_add = 0.2f64;
        let h2_gxe = 0.1f64;
        let env_sd = 1.5f64;
        let (a, b) = assign_effects(&mut rng, &p0, h2_add, h2_gxe, env_sd);
        let mut var_add = 0.0;
        let mut var_gxe = 0.0;
        for j in 0..l {
            let vj = 2.0 * p0[j] * (1.0 - p0[j]);
            var_add += a[j] * a[j] * vj;
            var_gxe += env_sd * env_sd * b[j] * b[j] * vj;
        }
        assert!((var_add - h2_add).abs() < 1e-9);
        assert!((var_gxe - h2_gxe).abs() < 1e-9);
    }

    #[test]
    fn test_compute_phenotypes_simple() {
        // p0 zero so dev = g
        let genotypes = vec![vec![0u8, 1u8, 2u8], vec![2u8, 0u8, 1u8]];
        let a = vec![1.0, -2.0, 0.5];
        let b = vec![0.0, 0.0, 0.0];
        let p0 = vec![0.0, 0.0, 0.0];
        let env = vec![0.0, 0.0];
        let z = compute_phenotypes(&genotypes, &a, &b, &p0, &env);
        // z0 = 0*1 + 1*(-2) + 2*0.5 = -2 + 1 = -1
        // z1 = 2*1 + 0*(-2) + 1*0.5 = 2 + 0.5 = 2.5
        assert!((z[0] + 1.0).abs() < 1e-9);
        assert!((z[1] - 2.5).abs() < 1e-9);
    }

    #[test]
    fn test_gaussian_fitness_normalization_and_order() {
        let phenos = vec![-1.0, 0.0, 1.0];
        let w = gaussian_fitness(&phenos, 0.0, 1.0);
        let sum: f64 = w.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert!(w[1] > w[0] && w[1] > w[2]);
        assert!((w[0] - w[2]).abs() < 1e-12);
    }

    #[test]
    fn test_allele_freqs_basic() {
        // 2 individuals, 2 loci
        let genotypes = vec![vec![0u8, 2u8], vec![1u8, 1u8]]; // totals: locus0=1, locus1=3
        let f = allele_freqs(&genotypes);
        assert!((f[0] - (1.0 / 4.0)).abs() < 1e-12);
        assert!((f[1] - (3.0 / 4.0)).abs() < 1e-12);
    }

    #[test]
    fn test_reintroduce_edges() {
        let mut rng = StdRng::seed_from_u64(3);
        let n = 5usize;
        // Case 1: all zero -> should introduce two copies
        let mut g1 = vec![vec![0u8; 1]; n];
        reintroduce_edges(&mut rng, &mut g1);
        let total1: i32 = g1.iter().map(|r| r[0] as i32).sum();
        assert_eq!(total1, 2);

        // Case 2: all two -> should reduce to 2N-2
        let mut g2 = vec![vec![2u8; 1]; n];
        reintroduce_edges(&mut rng, &mut g2);
        let total2: i32 = g2.iter().map(|r| r[0] as i32).sum();
        assert_eq!(total2, (2 * n) as i32 - 2);

        // Case 3: polymorphic -> unchanged total
        let mut g3 = vec![vec![0u8], vec![1u8], vec![2u8], vec![1u8], vec![0u8]];
        let before: i32 = g3.iter().map(|r| r[0] as i32).sum();
        reintroduce_edges(&mut rng, &mut g3);
        let after: i32 = g3.iter().map(|r| r[0] as i32).sum();
        assert_eq!(before, after);
    }

    #[test]
    fn test_reproduce_trivial_homozygous() {
        let mut rng = StdRng::seed_from_u64(4);
        // All parents 0 -> offspring all 0
        let parents0 = vec![vec![0u8; 3]; 10];
        let w = vec![1.0 / 10.0; 10];
        let kids0 = reproduce(&mut rng, &parents0, &w);
        assert!(kids0.iter().all(|row| row.iter().all(|&g| g == 0)));

        // All parents 2 -> offspring all 2
        let parents2 = vec![vec![2u8; 2]; 7];
        let w = vec![1.0 / 7.0; 7];
        let kids2 = reproduce(&mut rng, &parents2, &w);
        assert!(kids2.iter().all(|row| row.iter().all(|&g| g == 2)));
    }

    #[test]
    fn test_sample_index_basic() {
        let mut rng = StdRng::seed_from_u64(5);
        // Deterministic case: single-entry CDF always maps to index 0
        let cdf_single = vec![1.0];
        let uni = Uniform::new(0.0, 1.0);
        for _ in 0..10 {
            let idx = sample_index(&mut rng, &cdf_single, &uni);
            assert_eq!(idx, 0);
        }
    }

    #[test]
    fn test_draw_gamete_allele_values() {
        let mut rng = StdRng::seed_from_u64(6);
        assert_eq!(draw_gamete_allele(&mut rng, 0), 0);
        assert_eq!(draw_gamete_allele(&mut rng, 2), 1);
        for _ in 0..10 {
            let x = draw_gamete_allele(&mut rng, 1);
            assert!(x == 0 || x == 1);
        }
    }

    #[cfg(feature = "heavy-tests")]
    #[test]
    fn test_fixation_time_increases_with_population_size() {
        // In neutral drift, expected absorption time ~ O(N). Check monotonicity.
        let mut rng = StdRng::seed_from_u64(42);
        let reps = 100usize;
        let n_small = 50usize;
        let n_large = 200usize;
        let p0 = 0.5;
        let mut sum_small = 0u64;
        let mut sum_large = 0u64;
        for _ in 0..reps {
            let (t_s, _) = neutral_fix_or_loss_time(&mut rng, n_small, p0);
            sum_small += t_s;
            let (t_l, _) = neutral_fix_or_loss_time(&mut rng, n_large, p0);
            sum_large += t_l;
        }
        let avg_small = sum_small as f64 / reps as f64;
        let avg_large = sum_large as f64 / reps as f64;
        // Allow slack due to randomness; theoretical ratio ~ O(>~3x) for these Ns.
        assert!(
            avg_large > avg_small * 1.5,
            "expected avg_large ({avg_large}) > 1.5 * avg_small ({avg_small})"
        );
    }

    #[cfg(feature = "heavy-tests")]
    #[test]
    fn test_fixation_probability_matches_initial_frequency() {
        // Neutral fixation probability ≈ initial frequency p0
        let mut rng = StdRng::seed_from_u64(1337);
        let reps = 200usize;
        let n = 200usize;
        let p0 = 0.2;
        let mut fixed = 0usize;
        for _ in 0..reps {
            let (_, is_fixed) = neutral_fix_or_loss_time(&mut rng, n, p0);
            if is_fixed {
                fixed += 1;
            }
        }
        let p_fix = fixed as f64 / reps as f64;
        // Tolerance band to avoid flakiness in CI; 0.2 ± 0.1
        assert!(
            (p_fix - p0).abs() < 0.1,
            "fixation probability {p_fix:.3} not close to p0={p0}"
        );
    }
}
