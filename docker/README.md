# Containers

| Image | Dockerfile | What | Arch |
|---|---|---|---|
| `cookhla-rs` | `Dockerfile.cookhla-rs` | The **Rust port**, self-contained: `cookhla` + `beagle-rs` (imputation engine) + `plink` (from bioconda). Runs raw input → `.alleles`. | **amd64 + arm64** |
| `cookhla-oracle` | `Dockerfile.oracle` | The **original** CookHLA in its pinned conda env. The golden oracle — run it to produce reference output for parity. | amd64 |

The `cookhla-rs` image bundles everything it needs at runtime: the chr6 liftover chain is embedded
in the binary, `plink` comes from bioconda (published for both architectures), and `beagle-rs` is
built from the submodule. They are wired via `$PLINK` and `$BEAGLE_RS`.

## Run it

```sh
docker run --rm -v "$PWD/data:/data" ghcr.io/madhavajay/cookhla-rs \
    -i /data/1958BC.hg19 -hg 19 -o /data/out -ref /data/HM_CEU_REF \
    -gm /data/AGM.clpsB -ae /data/AGM.aver.erate
# -> /data/out.alleles
```

(`/data` is the working dir; mount your PLINK target, reference panel, and adaptive genetic map there.)

## Build it

Single arch (local, fast):

```sh
docker build -f docker/Dockerfile.cookhla-rs -t cookhla-rs .
```

Multi-arch (amd64 + arm64) — arm64 builds under QEMU emulation locally (slow); CI does the real
build/push:

```sh
docker/build-multiarch.sh                                   # build both, verify
PUSH=1 IMAGE=ghcr.io/madhavajay/cookhla-rs:latest docker/build-multiarch.sh   # build + push
PLATFORMS=linux/arm64 LOAD=1 docker/build-multiarch.sh      # one arch, load locally
```

CI (`.github/workflows/ci.yml`) builds the multi-arch image and pushes it to
`ghcr.io/<owner>/cookhla-rs` on every push to `main` and on tags.

## Generate the golden reference (oracle)

```sh
docker build -f docker/Dockerfile.oracle -t cookhla-oracle .
docker/oracle-run.sh          # -> fixtures/example/golden/...
docker/parity.sh              # diff cookhla-rs calls vs the golden
```

The example ships a precomputed adaptive genetic map, so MACH is not invoked for it — this
isolates the impute + consensus core. The `cookhla-rs` image produces HLA calls **identical** to
the oracle's on the example (verified on both amd64 and arm64).
