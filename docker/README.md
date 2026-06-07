# Containers

Two images, both built from the repo root so they can see `repos/CookHLA` (the reference
source, its vendored binaries, and the example data).

| Image | Dockerfile | What |
|---|---|---|
| `cookhla-oracle` | `Dockerfile.oracle` | The **original** CookHLA in its pinned conda env (`CookHLA_LINUX.yml`). The golden oracle — run it to produce reference outputs. |
| `cookhla-rs` | `Dockerfile.cookhla-rs` | The **Rust port** + the fast binaries it still shells out to (`plink`, `mach1`). Beagle is replaced in-process by `beagle-rs` (Phase 7). |

## Generate the golden reference

```sh
docker build -f docker/Dockerfile.oracle -t cookhla-oracle .
docker/oracle-run.sh          # -> fixtures/example/golden/1958BC+HM_CEU_REF.{alleles,hped}
```

The example ships a precomputed adaptive genetic map (`example/AGM.1958BC+HM_CEU_REF.*`), so
**MACH is not invoked** for this case — it isolates the impute + consensus core for first parity.

## Run the Rust port

```sh
docker build -f docker/Dockerfile.cookhla-rs -t cookhla-rs .
docker run --rm cookhla-rs --help
```

## Check parity

```sh
docker/parity.sh              # diffs HLA calls (must match) + probability drift (< 1e-4)
```

## Or via compose

```sh
docker compose -f docker/compose.yml build
docker compose -f docker/compose.yml run --rm oracle      # original CookHLA shell
docker compose -f docker/compose.yml run --rm cookhla-rs --help
```
