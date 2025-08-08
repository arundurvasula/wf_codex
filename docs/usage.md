# Usage & Install

## From source

Prerequisites:
- Rust toolchain (latest stable)

Commands:
- `cargo build` — compile
- `cargo run -- <args>` — run the binary
- `cargo test` — run tests
- `cargo fmt`, `cargo clippy -- -D warnings` — style & lints

## Prebuilt binaries

Download from GitHub Releases. Artifact names:
- `wf_codex-<platform>-<arch>` or `wf_codex-<platform>-<arch>-<cpu_tag>`
- Examples:
  - `wf_codex-linux-x86_64-x86-64-v3` — portable, fast on modern x86_64
  - `wf_codex-linux-x86_64` — native to CI runner CPU (may not run on older CPUs)
  - `wf_codex-macos-x86_64` — Intel Macs
  - `wf_codex-macos-arm64` — Apple Silicon Macs (M1/M2/M3)

Install example (macOS ARM64):
```bash
curl -L -o wf_codex "https://github.com/<org>/<repo>/releases/download/<tag>/wf_codex-macos-arm64"
chmod +x wf_codex
./wf_codex --help
```

Optional install to PATH:
```bash
sudo install -m 0755 wf_codex /usr/local/bin/wf_codex
```

## macOS Gatekeeper (unknown developer)

Since binaries are not notarized, macOS may block them. Options to allow:
- Finder: right-click → Open → Open
- System Settings → Privacy & Security → Security → “Open Anyway” (after first run attempt)
- Terminal: remove quarantine attribute
  - `xattr -d com.apple.quarantine ./wf_codex-macos-arm64`

## Quickstart

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

