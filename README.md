# CookHLA-rs

A fast Rust port of [CookHLA](https://github.com/WansonChoi/CookHLA) — accurate HLA imputation
(Cook et al., *Nat Commun* 2021). The default pipeline runs **raw genotypes → HLA calls in
~2 seconds** and reproduces the original CookHLA's calls **exactly** (verified against golden
output captured from the reference implementation).

```
cookhla -i target.hg19 -hg 19 -o OUT -ref PANEL -gm AGM.clpsB -ae AGM.aver.erate
# -> OUT.alleles  (FID IID gene 2-digit 4-digit pp1 pp2 conf)
```

## What it does

Reproduces CookHLA's default path entirely in compiled Rust, keeping the fast external binaries:

FixInput (liftover hg→18 + reference matching + de-ambiguation) → QC (`plink`) → reference panel
& adaptive-map building → CONVERT_IN → **9× Beagle imputation via [`beagle-rs`](https://github.com/madhavajay/beagle-rs)
(in parallel, no JVM)** → consensus caller → `.alleles`.

The interpreted layers of the original (R consensus, Python orchestration, Perl/csh/awk glue, and
the five Beagle utility JARs) are replaced by native Rust; `plink` is the only external tool kept
on the default path (`mach1`/MACH is only needed when the genetic map is auto-generated).

## Layout

| Path | What |
|---|---|
| `crates/cookhla` | library: file formats, converters, QC, panels, maps, consensus, pipeline |
| `crates/cookhla-cli` | the `cookhla` binary (argv-compatible with `CookHLA.py`) |
| `repos/CookHLA` | submodule: the original Python reference (+ vendored binaries + example data) |
| `repos/beagle-rs` | submodule: the Rust Beagle 5.5 port used as the imputation engine |
| `docker/` | reproducible oracle (original CookHLA) + cookhla-rs images |
| `TODO.md` | full porting status, per-stage parity results, remaining milestones |

## Build & test

```sh
git submodule update --init                              # CookHLA reference + beagle-rs
cargo build --release -p cookhla-cli                     # the `cookhla` binary
cargo build --release -p beagle-rs-cli \
    --manifest-path repos/beagle-rs/Cargo.toml           # the imputation engine
cargo test --workspace
```

Parity tests skip loudly unless the golden fixtures are present; regenerate them with
`docker build -f docker/Dockerfile.oracle -t cookhla-oracle . && docker/oracle-run.sh`.

## Docker

A self-contained image (`cookhla` + `beagle-rs` + `plink`) is built for **amd64 and arm64**:

```sh
docker run --rm -v "$PWD/data:/data" ghcr.io/madhavajay/cookhla-rs \
    -i /data/target.hg19 -hg 19 -o /data/out -ref /data/PANEL \
    -gm /data/AGM.clpsB -ae /data/AGM.aver.erate
```

See [`docker/README.md`](docker/README.md) and `docker/build-multiarch.sh`.

## Status

The default imputation path is complete and verified (raw input → identical HLA calls, 0
mismatches on the example). Not yet ported: MakeGeneticMap/MACH (only when `-gm`/`-ae` are
omitted), the legacy Beagle 4.1 / prephasing paths, and the accuracy module. See `TODO.md`.
