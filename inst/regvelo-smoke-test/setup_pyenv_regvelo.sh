#!/bin/bash
set -euo pipefail

# Create a dedicated pyenv + venv environment for the RegVelo GPU smoke test.
#
# This script is intended for clusters where software is managed with:
#   module load pyenv
# and pyenv is available, e.g. pyenv 2.5.2.
#
# It does not use conda or mamba.

PYTHON_VERSION="${PYTHON_VERSION:-3.10.14}"
VENV_DIR="${VENV_DIR:-${SCLET_SOURCE:-$PWD}/.venv-regvelo}"
PYENV_ROOT_USER="${PYENV_ROOT_USER:-$HOME/.pyenv-sclet-regvelo}"
TORCH_CUDA_INDEX="${TORCH_CUDA_INDEX:-https://download.pytorch.org/whl/cu121}"
PIP_INDEX_URL="${PIP_INDEX_URL:-}"
PIP_TRUSTED_HOST="${PIP_TRUSTED_HOST:-}"

module load pyenv

# Some clusters expose pyenv from a shared, read-only module path. python-build
# then tries to cache and install under that shared PYENV_ROOT and fails with
# permission errors. Force a user-writable pyenv root under $HOME.
export PYENV_ROOT="$PYENV_ROOT_USER"
mkdir -p "$PYENV_ROOT"
export PATH="$PYENV_ROOT/bin:$PYENV_ROOT/shims:$PATH"

if command -v pyenv >/dev/null 2>&1; then
  eval "$(pyenv init -)"
else
  echo "pyenv command not found after module load pyenv" >&2
  exit 1
fi

echo "Using PYENV_ROOT: $PYENV_ROOT"

if ! pyenv versions --bare | grep -qx "$PYTHON_VERSION"; then
  echo "Installing Python $PYTHON_VERSION with pyenv..."
  pyenv install "$PYTHON_VERSION"
else
  echo "Python $PYTHON_VERSION already exists in pyenv."
fi

PYTHON_BIN="$(pyenv root)/versions/$PYTHON_VERSION/bin/python"
echo "Using Python: $PYTHON_BIN"

"$PYTHON_BIN" -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

python -m pip install --upgrade pip setuptools wheel

PIP_COMMON_ARGS=()
if [[ -n "$PIP_INDEX_URL" ]]; then
  PIP_COMMON_ARGS+=(--index-url "$PIP_INDEX_URL")
fi
if [[ -n "$PIP_TRUSTED_HOST" ]]; then
  PIP_COMMON_ARGS+=(--trusted-host "$PIP_TRUSTED_HOST")
fi

# Install heavy compiled dependencies from prebuilt wheels only. If a wheel is
# unavailable, fail early rather than falling back to source builds that may
# require Bazel, compiler toolchains, or external downloads with broken CA
# certificates.
python -m pip install ${PIP_COMMON_ARGS[@]+"${PIP_COMMON_ARGS[@]}"} --only-binary=:all: \
  "numpy<2" \
  "scipy" \
  "pandas>=2.2" \
  "h5py>=3.10" \
  "tensorstore>=0.1.45" \
  "ml-dtypes" \
  "jax>=0.4.20" \
  "jaxlib>=0.4.20"

# RegVelo metadata expects torch < 2.6. Use a CUDA-enabled 2.5.x wheel.
python -m pip install \
  --index-url "$TORCH_CUDA_INDEX" \
  "torch==2.5.1" \
  "torchvision==0.20.1" \
  "torchaudio==2.5.1"

python -m pip install ${PIP_COMMON_ARGS[@]+"${PIP_COMMON_ARGS[@]}"} --prefer-binary \
  "anndata>=0.10.8" \
  "scanpy>=1.10.3" \
  "scvelo>=0.3.2" \
  "cellrank==2.0.7" \
  "scvi-tools>=1.0,<1.2.1" \
  "regvelo==0.4.2"

python - <<'PY'
import sys
import torch
import regvelo
import scvelo
import cellrank
import scvi
import scanpy
import anndata

print("python", sys.executable)
print("torch", torch.__version__)
print("cuda_available", torch.cuda.is_available())
print("cuda_device_count", torch.cuda.device_count())
print("regvelo", getattr(regvelo, "__version__", "unknown"))
print("scvelo", scvelo.__version__)
print("cellrank", cellrank.__version__)
print("scvi", scvi.__version__)
print("scanpy", scanpy.__version__)
print("anndata", anndata.__version__)
PY

echo "Created environment: $VENV_DIR"
echo "Activate with:"
echo "  source $VENV_DIR/bin/activate"
