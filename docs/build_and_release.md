# Build & Release

## Local build

```bash
cargo build --release
```

Release profile is tuned for runtime speed:
- fat LTO, 1 codegen unit, opt-level=3, panic=abort, strip symbols

## CI builds & artifacts

The GitHub Actions release workflow produces:
- Linux (native): `wf_codex-linux-x86_64` (tuned to runner CPU)
- Linux (portable, fast): `wf_codex-linux-x86_64-x86-64-v3`
- macOS Intel: `wf_codex-macos-x86_64`
- macOS Apple Silicon: `wf_codex-macos-arm64`

For testing (manual workflow_dispatch), artifacts are uploaded to the workflow run.
For published releases, assets are attached to the GitHub Release.

## API docs

Generate locally:
```bash
cargo doc --no-deps --open
```

If enabled, CI publishes `cargo doc` output to GitHub Pages.

