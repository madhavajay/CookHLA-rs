#!/usr/bin/env bash
# Run the ORIGINAL CookHLA on the bundled example inside the oracle container and copy the
# golden outputs out. This is the parity reference cookhla-rs must reproduce.
#
# Usage (from repo root, after `docker build -f docker/Dockerfile.oracle -t cookhla-oracle .`):
#   docker/oracle-run.sh                # -> fixtures/example/golden/
#
# The example ships a precomputed AGM (example/AGM.1958BC+HM_CEU_REF.*), so MACH is NOT
# invoked here — this isolates the core impute+consensus pipeline for first parity.
set -euo pipefail

IMAGE="${IMAGE:-cookhla-oracle}"
OUTDIR="${OUTDIR:-$(cd "$(dirname "$0")/.." && pwd)/fixtures/example/golden}"
mkdir -p "$OUTDIR"

# The reference command from repos/CookHLA/README.md, writing into a container tmp dir; we then
# copy the *.alleles / *.hped / accuracy back to the host golden dir.
docker run --rm -v "$OUTDIR:/golden" "$IMAGE" -c '
  set -euo pipefail
  cd /opt/CookHLA
  mkdir -p /tmp/run
  python CookHLA.py \
      -i example/1958BC.hg19 \
      -hg 19 \
      -o /tmp/run/1958BC+HM_CEU_REF \
      -ref example/HM_CEU_REF \
      -gm example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB \
      -ae example/AGM.1958BC+HM_CEU_REF.aver.erate \
      -mem 2g
  echo "--- outputs ---"; ls -la /tmp/run
  cp /tmp/run/1958BC+HM_CEU_REF.alleles /golden/ 2>/dev/null || true
  cp /tmp/run/1958BC+HM_CEU_REF.hped    /golden/ 2>/dev/null || true
'
echo "Golden outputs in: $OUTDIR"
ls -la "$OUTDIR"
