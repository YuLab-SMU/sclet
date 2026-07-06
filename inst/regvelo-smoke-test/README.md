# sclet RegVelo GPU Smoke Test

This directory is a self-contained server-side smoke-test kit for the RegVelo
integration in `sclet`. It is separate from package runtime code and from the
bookdown manuscript.

The test verifies that a user-writable pyenv + venv environment can import the
required Python stack, CUDA is visible to PyTorch inside SLURM, `RunRegVelo()`
writes a velocity assay into a `SingleCellExperiment`, and `RunCellRank()` can
consume that stored RegVelo velocity assay.

## Files

- `setup_pyenv_regvelo.sh`: creates a private pyenv + venv Python environment
  without conda/mamba.
- `download_smoke_h5ad.py`: downloads the public scVelo pancreas h5ad once.
- `prepare_data.R`: reads a local h5ad through Python `anndata` and writes a
  small SCE subset plus a synthetic prior GRN.
- `run_regvelo_smoke.R`: runs `RunRegVelo()` and writes the first-stage RDS.
- `run_cellrank_after_regvelo.R`: reads the first-stage RDS, runs CellRank, and
  writes diagnostics/plots.
- `analyze_regvelo_result.R`: generates quick RegVelo velocity magnitude
  summaries from the first-stage RDS.
- `slurm_download_h5ad.yulab-pyenv.sbatch`: YuLab-style h5ad download job.
- `slurm_prepare_data.yulab-pyenv.sbatch`: YuLab-style data-preparation job.
- `slurm_run_regvelo.yulab-gpu-pyenv.sbatch`: YuLab-style GPU RegVelo job.
- `slurm_run_cellrank_after_regvelo.yulab-gpu-pyenv.sbatch`: YuLab-style
  CellRank-after-RegVelo job.

## 1. Configure paths

```bash
export SCLET_SOURCE=/path/to/sclet
cd "$SCLET_SOURCE"
```

The R scripts set `~/Rlib` first while preserving existing `.libPaths()`.

## 2. Create the Python environment

```bash
cd "$SCLET_SOURCE/inst/regvelo-smoke-test"
PYENV_ROOT_USER="$HOME/.pyenv-sclet-regvelo" \
VENV_DIR="$SCLET_SOURCE/.venv-regvelo" \
bash setup_pyenv_regvelo.sh
```

This installs Python 3.10.14 under `$HOME/.pyenv-sclet-regvelo` and creates
`$SCLET_SOURCE/.venv-regvelo`. It installs PyTorch 2.5.1 CUDA wheels plus the
Python packages needed by RegVelo and CellRank.

If the cluster uses an internal PyPI mirror or has SSL interception:

```bash
PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple \
PIP_TRUSTED_HOST=pypi.tuna.tsinghua.edu.cn \
PYENV_ROOT_USER="$HOME/.pyenv-sclet-regvelo" \
VENV_DIR="$SCLET_SOURCE/.venv-regvelo" \
bash "$SCLET_SOURCE/inst/regvelo-smoke-test/setup_pyenv_regvelo.sh"
```

## 3. Prepare the source h5ad

If the server can download the public scVelo dataset:

```bash
cd "$SCLET_SOURCE"
mkdir -p logs
sbatch inst/regvelo-smoke-test/slurm_download_h5ad.yulab-pyenv.sbatch
```

If GitHub access is unreliable, download this file locally and upload it to
`regvelo_smoke_source/pancreas.h5ad`:

```text
https://raw.githubusercontent.com/theislab/scvelo_notebooks/master/data/Pancreas/endocrinogenesis_day15.h5ad
```

## 4. Build the small smoke-test dataset

Interactive R route:

```r
source("inst/regvelo-smoke-test/prepare_data.R")
prepare_regvelo_smoke_data(
  source_h5ad = "regvelo_smoke_source/pancreas.h5ad",
  outdir = "regvelo_smoke_data",
  n_cells = 300,
  n_genes = 800
)
```

SLURM route:

```bash
cd "$SCLET_SOURCE"
mkdir -p logs
sbatch inst/regvelo-smoke-test/slurm_prepare_data.yulab-pyenv.sbatch
```

Expected outputs:

- `regvelo_smoke_data/sce_regvelo_smoke.rds`
- `regvelo_smoke_data/prior_grn.csv`
- `regvelo_smoke_data/manifest.csv`

## 5. Run RegVelo on GPU

```bash
cd "$SCLET_SOURCE"
mkdir -p logs
PYENV_ROOT_USER="$HOME/.pyenv-sclet-regvelo" \
VENV_DIR="$SCLET_SOURCE/.venv-regvelo" \
sbatch inst/regvelo-smoke-test/slurm_run_regvelo.yulab-gpu-pyenv.sbatch
```

Expected outputs:

- `regvelo_smoke_results/sce_regvelo_smoke_result.rds`
- `regvelo_smoke_results/summary.csv`
- `logs/regvelo-yulab-<jobid>.out`
- `logs/regvelo-yulab-<jobid>.err`

The summary should show `torch_cuda_available = TRUE` and a velocity assay such
as `regvelo_smoke_velocity`.

## 6. Run CellRank from the RegVelo velocity assay

This second-stage job reads the RegVelo result RDS and does not retrain RegVelo.

```bash
cd "$SCLET_SOURCE"
mkdir -p logs
PYENV_ROOT_USER="$HOME/.pyenv-sclet-regvelo" \
VENV_DIR="$SCLET_SOURCE/.venv-regvelo" \
sbatch inst/regvelo-smoke-test/slurm_run_cellrank_after_regvelo.yulab-gpu-pyenv.sbatch
```

The smoke result does not contain a biological cell-type column. If
`CLUSTER_KEY` is not supplied, the script creates a temporary
`regvelo_smoke_cluster` column using k-means on PCA. For real data, pass a real
annotation:

```bash
CLUSTER_KEY=celltype \
sbatch inst/regvelo-smoke-test/slurm_run_cellrank_after_regvelo.yulab-gpu-pyenv.sbatch
```

Expected outputs:

- `regvelo_smoke_cellrank/sce_regvelo_cellrank_result.rds`
- `regvelo_smoke_cellrank/summary.csv`
- `regvelo_smoke_cellrank/cellrank_terminal_state_counts.csv`
- `regvelo_smoke_cellrank/cellrank_fate_probability_summary.csv`
- `regvelo_smoke_cellrank/velocity_fate_correlation.csv`
- `regvelo_smoke_cellrank/cellrank_terminal_states.png`
- `regvelo_smoke_cellrank/cellrank_fate_probability_1.png`
- `regvelo_smoke_cellrank/velocity_fate_correlation.png`
- `logs/cellrank-yulab-<jobid>.out`
- `logs/cellrank-yulab-<jobid>.err`

For the small smoke-test dataset, a single terminal state or near-constant fate
probability is acceptable. It means the pipeline is connected; it should not be
treated as biological evidence for a real branching process.

## 7. R module troubleshooting

The YuLab SLURM templates try `R`, `R/4.4.0`, and `R/4.3.0` as module names. If
your cluster uses a different name, pass it explicitly:

```bash
R_MODULE_CANDIDATES="R/4.4.0 R" \
sbatch inst/regvelo-smoke-test/slurm_run_cellrank_after_regvelo.yulab-gpu-pyenv.sbatch
```

If module loading is unreliable but the executable path is known:

```bash
R_BIN=/biostack/tools/devtools/R/4.4.0/bin/R \
sbatch inst/regvelo-smoke-test/slurm_run_cellrank_after_regvelo.yulab-gpu-pyenv.sbatch
```

The templates avoid `module purge` by default because cluster R builds often
need compiler runtime libraries already present in the module stack.

## 8. What to send back if it fails

Send the relevant command plus the matching SLURM logs and summary files:

- `logs/regvelo-yulab-<jobid>.out`
- `logs/regvelo-yulab-<jobid>.err`
- `logs/cellrank-yulab-<jobid>.out`
- `logs/cellrank-yulab-<jobid>.err`
- `regvelo_smoke_results/summary.csv`, if present
- `regvelo_smoke_cellrank/summary.csv`, if present
