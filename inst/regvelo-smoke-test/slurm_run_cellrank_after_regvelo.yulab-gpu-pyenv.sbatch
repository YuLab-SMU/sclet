#!/bin/bash
#SBATCH -p gpu_batch
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 04:00:00
#SBATCH --mem-per-cpu 8000
#SBATCH --job-name sclet-cellrank
#SBATCH --output logs/cellrank-yulab-%J.out
#SBATCH --error logs/cellrank-yulab-%J.err

set -euo pipefail

# Site-specific helper for the YuLab/SMU module + pyenv server. This is a
# second-stage smoke test: it consumes the RDS produced by
# run_regvelo_smoke.R and does not retrain RegVelo.

SCLET_SOURCE="${SCLET_SOURCE:-$PWD}"
SMOKE_DIR="${SMOKE_DIR:-$SCLET_SOURCE/inst/regvelo-smoke-test}"
INPUT_RDS="${INPUT_RDS:-$PWD/regvelo_smoke_results/sce_regvelo_smoke_result.rds}"
OUT_DIR="${OUT_DIR:-$PWD/regvelo_smoke_cellrank}"
VENV_DIR="${VENV_DIR:-$SCLET_SOURCE/.venv-regvelo}"
PYENV_ROOT_USER="${PYENV_ROOT_USER:-$HOME/.pyenv-sclet-regvelo}"
VELOCITY_ID="${VELOCITY_ID:-regvelo_smoke}"
TRAJECTORY_ID="${TRAJECTORY_ID:-cellrank_regvelo_smoke}"
REDUCTION="${REDUCTION:-PCA}"
CLUSTER_KEY="${CLUSTER_KEY:-}"
K="${K:-4}"
PURGE_MODULES="${PURGE_MODULES:-0}"
R_MODULE="${R_MODULE:-R}"
R_MODULE_CANDIDATES="${R_MODULE_CANDIDATES:-$R_MODULE R/4.4.0 R/4.3.0}"
PYENV_MODULE="${PYENV_MODULE:-pyenv}"
EXTRA_MODULES="${EXTRA_MODULES:-}"
R_BIN="${R_BIN:-}"
RSCRIPT_BIN="${RSCRIPT_BIN:-}"

mkdir -p logs "$OUT_DIR"

if ! type module >/dev/null 2>&1; then
  for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
    if [[ -r "$init" ]]; then
      # shellcheck disable=SC1090
      source "$init"
      break
    fi
  done
fi

if [[ "$PURGE_MODULES" == "1" ]]; then
  module purge
fi
for mod in $EXTRA_MODULES; do
  module load "$mod"
done
module load "$PYENV_MODULE"

load_r_module() {
  if [[ -n "$R_BIN" && -x "$R_BIN" ]]; then
    return 0
  fi
  if command -v R >/dev/null 2>&1 || command -v Rscript >/dev/null 2>&1; then
    return 0
  fi
  local mod
  for mod in $R_MODULE_CANDIDATES; do
    module load "$mod" >/dev/null 2>&1 || continue
    if command -v R >/dev/null 2>&1 || command -v Rscript >/dev/null 2>&1; then
      echo "Loaded R module: $mod"
      return 0
    fi
  done
  echo "Could not find R after trying R_MODULE_CANDIDATES='$R_MODULE_CANDIDATES'." >&2
  echo "Loaded modules:" >&2
  module list >&2 || true
  echo "Available R-like modules:" >&2
  module avail R >&2 || true
  return 1
}
load_r_module
export PYENV_ROOT="$PYENV_ROOT_USER"
export PATH="$PYENV_ROOT/bin:$PYENV_ROOT/shims:$PATH"
eval "$(pyenv init -)"
source "$VENV_DIR/bin/activate"

run_r_script() {
  local script="$1"
  shift
  if [[ -n "$RSCRIPT_BIN" && -x "$RSCRIPT_BIN" ]]; then
    "$RSCRIPT_BIN" "$script" "$@"
  elif command -v Rscript >/dev/null 2>&1; then
    "$(command -v Rscript)" "$script" "$@"
  elif [[ -n "$R_BIN" && -x "$R_BIN" ]]; then
    "$R_BIN" --vanilla --slave -f "$script" --args "$@"
  elif command -v R >/dev/null 2>&1; then
    "$(command -v R)" --vanilla --slave -f "$script" --args "$@"
  else
    echo "Neither Rscript nor R was found after module loading." >&2
    echo "Set R_MODULE, R_MODULE_CANDIDATES, R_BIN, or RSCRIPT_BIN when submitting." >&2
    exit 1
  fi
}

echo "Running on host: $(hostname)"
module list || true
echo "Python: $(which python)"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-unset}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi || true

echo "R: $(command -v R || true)"
echo "Rscript: $(command -v Rscript || true)"

run_r_script "$SMOKE_DIR/run_cellrank_after_regvelo.R" \
  --sclet-path="$SCLET_SOURCE" \
  --result="$INPUT_RDS" \
  --outdir="$OUT_DIR" \
  --python="$(which python)" \
  --velocity-id="$VELOCITY_ID" \
  --trajectory-id="$TRAJECTORY_ID" \
  --reduction="$REDUCTION" \
  --cluster-key="$CLUSTER_KEY" \
  --k="$K"
