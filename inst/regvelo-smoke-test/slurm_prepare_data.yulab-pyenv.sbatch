#!/bin/bash
#SBATCH -p gpu_batch
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH --mem-per-cpu 8000
#SBATCH --job-name sclet-regvelo-data
#SBATCH --output logs/regvelo-data-yulab-%J.out
#SBATCH --error logs/regvelo-data-yulab-%J.err

set -euo pipefail

# Site-specific template for servers that use `module load pyenv`, matching the
# provided Jupyter GPU job style. This template uses a plain Python venv created
# by setup_pyenv_regvelo.sh, so it does not require pyenv-virtualenv.

SCLET_SOURCE="${SCLET_SOURCE:-$PWD}"
SMOKE_DIR="${SMOKE_DIR:-$SCLET_SOURCE/inst/regvelo-smoke-test}"
DATA_DIR="${DATA_DIR:-$PWD/regvelo_smoke_data}"
VENV_DIR="${VENV_DIR:-$SCLET_SOURCE/.venv-regvelo}"
PYENV_ROOT_USER="${PYENV_ROOT_USER:-$HOME/.pyenv-sclet-regvelo}"
SOURCE_H5AD="${SOURCE_H5AD:-}"
PURGE_MODULES="${PURGE_MODULES:-0}"
R_MODULE="${R_MODULE:-R}"
PYENV_MODULE="${PYENV_MODULE:-pyenv}"
EXTRA_MODULES="${EXTRA_MODULES:-}"

mkdir -p logs "$DATA_DIR"

if [[ "$PURGE_MODULES" == "1" ]]; then
  module purge
fi
for mod in $EXTRA_MODULES; do
  module load "$mod"
done
module load "$PYENV_MODULE"
module load "$R_MODULE"
export PYENV_ROOT="$PYENV_ROOT_USER"
export PATH="$PYENV_ROOT/bin:$PYENV_ROOT/shims:$PATH"
eval "$(pyenv init -)"
source "$VENV_DIR/bin/activate"

echo "Running on host: $(hostname)"
module list || true
echo "Python: $(which python)"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-unset}"
echo "R: $(command -v R || true)"
echo "Rscript: $(command -v Rscript || true)"

if [[ -n "$SOURCE_H5AD" ]]; then
  R --vanilla -e "source(file.path('$SMOKE_DIR', 'prepare_data.R')); prepare_regvelo_smoke_data(source_h5ad = '$SOURCE_H5AD', outdir = '$DATA_DIR', n_cells = 300, n_genes = 800)"
else
  R --vanilla -e "source(file.path('$SMOKE_DIR', 'prepare_data.R')); prepare_regvelo_smoke_data(outdir = '$DATA_DIR', n_cells = 300, n_genes = 800)"
fi
