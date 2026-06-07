#!/usr/bin/env bash
# Run cookhla-rs on the example and diff its HLA calls against the oracle's golden output.
#
# Parity bar (v1): identical allele *calls* in the *.alleles file; posterior-probability columns
# may differ within a small epsilon (MACH RNG + FP accumulation). This script does the exact
# call comparison and reports probability drift. Wire-up completes with Phase 7 (impute) and
# Phase 5 (consensus); until then it documents the intended gate.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GOLDEN="${GOLDEN:-$ROOT/fixtures/example/golden/1958BC+HM_CEU_REF.alleles}"
RS_OUT="${RS_OUT:-$ROOT/fixtures/example/rs/1958BC+HM_CEU_REF.alleles}"

if [[ ! -f "$GOLDEN" ]]; then
  echo "ERROR: golden not found ($GOLDEN). Run docker/oracle-run.sh first." >&2
  exit 1
fi
if [[ ! -f "$RS_OUT" ]]; then
  echo "ERROR: cookhla-rs output not found ($RS_OUT). Run the Rust pipeline first." >&2
  exit 1
fi

# Columns of *.alleles: FID IID gene 1-digit(allele1,allele2) 4-digit(a1,a2) pp1 pp2 conf.
# Calls = cols 1-5; probabilities = cols 6-8.
calls() { awk '{print $1,$2,$3,$4,$5}' "$1"; }

echo "== HLA call parity (cols 1-5, must be identical) =="
if diff <(calls "$GOLDEN") <(calls "$RS_OUT"); then
  echo "PASS: HLA calls identical"
else
  echo "FAIL: HLA calls differ" >&2
  exit 1
fi

echo "== probability drift (cols 6-8, must be < 1e-4) =="
paste "$GOLDEN" "$RS_OUT" | awk '
  { for (k=6;k<=8;k++) { d=($k)-($(k+8)); if (d<0) d=-d; if (d>max) max=d } }
  END { printf "max |Δprob| = %.3e\n", max; if (max > 1e-4) { print "FAIL: drift too large" > "/dev/stderr"; exit 1 } print "PASS" }'
