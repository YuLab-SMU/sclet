#!/bin/bash
#SBATCH -p gpu_batch
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH --mem-per-cpu 8000
#SBATCH --job-name sclet-regvelo-h5ad
#SBATCH --output logs/regvelo-h5ad-yulab-%J.out
#SBATCH --error logs/regvelo-h5ad-yulab-%J.err

set -euo pipefail

# Site-specific template for servers that use `module load pyenv`, matching the
# provided Jupyter GPU job style. This step only downloads/caches h5ad data and
# does not require a GPU, but it uses the same pyenv/venv environment.

SCLET_SOURCE="${SCLET_SOURCE:-$PWD}"
SMOKE_DIR="${SMOKE_DIR:-$SCLET_SOURCE/inst/regvelo-smoke-test}"
VENV_DIR="${VENV_DIR:-$SCLET_SOURCE/.venv-regvelo}"
PYENV_ROOT_USER="${PYENV_ROOT_USER:-$HOME/.pyenv-sclet-regvelo}"
SOURCE_H5AD="${SOURCE_H5AD:-$SCLET_SOURCE/regvelo_smoke_source/pancreas.h5ad}"
PURGE_MODULES="${PURGE_MODULES:-0}"
PYENV_MODULE="${PYENV_MODULE:-pyenv}"
EXTRA_MODULES="${EXTRA_MODULES:-}"

mkdir -p logs "$(dirname "$SOURCE_H5AD")"

if [[ "$PURGE_MODULES" == "1" ]]; then
  module purge
fi
for mod in $EXTRA_MODULES; do
  module load "$mod"
done
module load "$PYENV_MODULE"
export PYENV_ROOT="$PYENV_ROOT_USER"
export PATH="$PYENV_ROOT/bin:$PYENV_ROOT/shims:$PATH"
eval "$(pyenv init -)"
source "$VENV_DIR/bin/activate"

echo "Running on host: $(hostname)"
module list || true
echo "Python: $(which python)"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-unset}"
echo "Output h5ad: $SOURCE_H5AD"

python "$SMOKE_DIR/download_smoke_h5ad.py" --out "$SOURCE_H5AD"
