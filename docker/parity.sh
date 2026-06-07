#!/usr/bin/env bash
# Run cookhla-rs on the bundled example and diff its HLA calls against the oracle golden.
#
# Parity bar: the allele *calls* (cols 1-5 of *.alleles) must be identical; the posterior columns
# may differ slightly (beagle-rs is Beagle 5.5 vs the oracle's 5.1).
#
#   docker/parity.sh                       # uses the local `cookhla-rs` image
#   IMAGE=ghcr.io/madhavajay/cookhla-rs docker/parity.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE="${IMAGE:-cookhla-rs}"
EXAMPLE="$ROOT/repos/CookHLA/example"
GOLDEN="$ROOT/fixtures/example/golden/1958BC+HM_CEU_REF.MHC.HLA_IMPUTATION_OUT.alleles"

[[ -f "$GOLDEN" ]] || { echo "golden not found ($GOLDEN). Run docker/oracle-run.sh first." >&2; exit 1; }

out="$(mktemp -d)"
docker run --rm -v "$EXAMPLE:/data:ro" -v "$out:/out" "$IMAGE" \
    -i /data/1958BC.hg19 -hg 19 -o /out/result -ref /data/HM_CEU_REF \
    -gm /data/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB \
    -ae /data/AGM.1958BC+HM_CEU_REF.aver.erate

echo "== HLA call parity (FID IID gene 2-digit 4-digit) =="
if diff <(awk '{print $1,$2,$3,$4,$5}' "$out/result.alleles") \
        <(awk '{print $1,$2,$3,$4,$5}' "$GOLDEN"); then
  echo "PASS: HLA calls identical to the oracle golden"
else
  echo "FAIL: HLA calls differ" >&2
  exit 1
fi
