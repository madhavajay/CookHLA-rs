#!/usr/bin/env bash
# Build (and optionally push) the multi-arch cookhla-rs image (linux/amd64 + linux/arm64).
#
#   docker/build-multiarch.sh                         # build both arches (verify only, no load)
#   PUSH=1 IMAGE=ghcr.io/madhavajay/cookhla-rs:latest docker/build-multiarch.sh   # build + push
#   PLATFORMS=linux/arm64 LOAD=1 docker/build-multiarch.sh   # build one arch and load locally
#
# arm64 on an amd64 host builds under QEMU emulation — slow (the Rust compile is emulated). CI
# (.github/workflows/ci.yml) does the real multi-arch build/push to ghcr.io.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE="${IMAGE:-cookhla-rs}"
PLATFORMS="${PLATFORMS:-linux/amd64,linux/arm64}"

# Register cross-arch emulation + a container builder (idempotent).
docker run --privileged --rm tonistiigi/binfmt --install arm64,amd64 >/dev/null 2>&1 || true
docker buildx inspect cookhla-builder >/dev/null 2>&1 \
  || docker buildx create --name cookhla-builder --driver docker-container --bootstrap >/dev/null

if [[ "${PUSH:-0}" == "1" ]]; then
  out=(--push)
elif [[ "${LOAD:-0}" == "1" ]]; then
  out=(--load) # single-platform only
else
  out=(--output type=cacheonly) # build to verify, don't load a manifest list locally
fi

set -x
docker buildx build --builder cookhla-builder \
  --platform "$PLATFORMS" \
  -f "$ROOT/docker/Dockerfile.cookhla-rs" \
  -t "$IMAGE" "${out[@]}" "$ROOT"
