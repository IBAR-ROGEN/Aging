#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mkdir -p .tools
ARCH="$(uname -m)"
case "${ARCH}" in
arm64) PLATFORM="osx-arm64" ;;
x86_64) PLATFORM="osx-64" ;;
*)
  echo "Unsupported CPU architecture: ${ARCH}" >&2
  exit 1
  ;;
esac
if [ ! -x .tools/bin/micromamba ]; then
  curl -fsSL "https://micro.mamba.pm/api/micromamba/${PLATFORM}/latest" | tar -xj -C .tools bin/micromamba
  chmod +x .tools/bin/micromamba
fi
.tools/bin/micromamba create -y -p "${ROOT}/.r-env" -c conda-forge r-base
