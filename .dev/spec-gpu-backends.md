# Spec: GPU acceleration for the torch-backed Python methods

> Status: **draft for review** (2026-10-23)
> Baseline: `devel`, package 1.1.0, uncommitted work in the tree.
> Decision taken before this spec: **the default stays CPU. GPU is opt-in.**
> Rationale: an environment carrying CUDA wheels is ~3 GB instead of ~200 MB per
> environment, and a package whose default install downloads 15 GB of CUDA
> libraries to serve a minority of users is a package nobody installs.

## 1. Problem

Five of sclet's seven Python backends are (or wrap) deep-learning models that
run on a GPU one or two orders of magnitude faster than on a CPU:

| env | backend | torch? | measured benefit |
|---|---|---|---|
| `sclet_scvi_env` | scvi-tools 1.4.3 (`RunIntegration(method="scVI")`) | yes (via pytorch-lightning) | the biggest one on this list |
| `sclet_cell2location_env` | cell2location 0.1.5 (`RunSpatialDeconvolution`) | yes (via scvi-tools) | large |
| `sclet_regvelo_env` | regvelo 0.4.2 (`RunRegVelo`) | yes (`pytorch=2.5.1`, already pinned) | large |
| `sclet_celloracle_env` | celloracle 0.20.0 (`RunInSilicoPerturbation`) | no — CellOracle is numpy/sklearn | **none** |
| `sclet_cellrank_env` | cellrank 2.0.7 + scvelo 0.3.4 (`RunCellRank`, `RunVelocity`) | no (CPU eigendecomposition) | none, unless scVelo's `recover_dynamics` is used |
| `sclet_scenic_env` | pySCENIC (`RunSCENIC`) | no (arboreto boosting, dask-parallel) | CPU + threads |
| `sclet_cellphonedb_env` | CellPhoneDB (`RunCellPhoneDB`) | no | none |

Verified on this machine (RTX 4090, torch 2.13.0+cu130): a 4000³ matmul ×20 takes
**0.28 s on GPU vs 2.21 s on CPU**. And in `~/data/dl-learning/.venv` a
`scvi-tools + cellrank + scvelo` resolution **reuses the venv's CUDA torch**
(`torch` untouched, only `torchmetrics` added), which is the evidence that a
pip/venv basilisk environment can carry GPU support without conda.

Today none of this is reachable: `CUDA_VISIBLE_DEVICES`, `accelerator` and
`device` appear **nowhere** in `R/`, and `R/regvelo.R:24-25` tells the user that
GPU requires `backend = "reticulate"` and their own Python — i.e. the package
pushes its own GPU support onto the user.

## 2. Design

### 2.1 One new environment, not seven

`sclet_gpu_env` — a single GPU-capable environment holding the *shared* GPU
stack (CUDA torch + the scvi family), used by the methods that benefit. It is
declared alongside the existing seven and keyed on its own spec hash by the
existing `sclet_env_spec_fingerprint()` machinery, so adding it does not move any
existing environment path.

What it pins (all values chosen against what this machine already resolves):

```
python=3.12
scanpy=1.12.3      # matches ~/data/dl-learning
anndata=0.13.2
numpy=2.4.4
pandas=3.0.5
torch=2.13.0
pytorch-cuda=13.0  # the wheel channel that produced torch 2.13.0+cu130
scvi-tools=1.5.1
scvelo=0.3.4
cellrank=2.3.3
```

Not in it, with reasons: **celloracle** (no torch), **pySCENIC/CellPhoneDB**
(no torch), **regvelo** (its own env keeps `pytorch=2.5.1` because regvelo 0.4.2
pins a different scvi-tools; a GPU variant of regvelo is §8 Q2).

### 2.2 How a method picks an environment

One resolver, used by the three methods that gain a GPU path:

```r
sclet_ai_device <- function(device = c("auto", "cpu", "cuda")) { ... }
sclet_torch_backend <- function(device, method_fn, gpu_fn) { ... }
```

- `device = "cpu"` (default) → the existing environment and the existing code
  path. **No behaviour change for anyone who does not ask.**
- `device = "cuda"` → `sclet_gpu_env` and a GPU code path.
- `device = "auto"` → GPU if the environment exists and reports a device,
  else CPU with one `message()` saying so.

The methods gain a `device` parameter: `RunIntegration()`, `RunSpatialDeconvolution()`
(and the `run_scvi`-style action later). `RunCellRank`/`RunVelocity` gain nothing
in this milestone (no torch), which is stated rather than left ambiguous.

### 2.3 Inside the GPU callbacks

Each callback takes the device from its arguments and applies it the way its
library expects, never unconditionally:

| backend | how |
|---|---|
| scvi-tools | `SCVI.setup_anndata(..., accelerator=..., devices=...)` or `model.train(accelerator=...)`; scvi-tools 1.5 accepts `accelerator="gpu"`/`"cpu"`/`"auto"` |
| cell2location | passes its scvi-tools model the same way |
| cellrank / scvelo | not touched (CPU) |

A device that is requested but unavailable raises a **typed error**
(`sclet_gpu_unavailable`) naming `sclet_gpu_status()` — not a Python traceback
from inside `torch.cuda`.

### 2.4 First-class GPU introspection

```r
sclet_gpu_status()   # data.frame: available, torch, cuda, device names, env path
```

Reads the `sclet_gpu_env` interpreter (building nothing) and reports what it
sees. This is the "is my card even visible to the environment sclet would use"
question, answerable before a multi-hour run.

### 2.5 How a user turns it on

Three levels, all opt-in:

```r
RunIntegration(sce, method = "scVI", batch = "batch", device = "cuda")   # per call
options(sclet.python.device = "cuda")                                   # per session
Sys.setenv(SCLET_PYTHON_DEVICE = "cuda")                                # per machine
```

The option wins over the default, the argument wins over the option, and the
environment variable is the fallback for non-interactive use. `"auto"` is
accepted everywhere.

### 2.6 Environment construction cost

`sclet_gpu_env` is built lazily, like the other seven: the first
`device = "cuda"` call pays ~3 GB and several minutes, and `render`/NEWS/book
say so. A user who never asks never pays. If the build fails, the existing
`sclet_env_setup_error()` path already carries repair text.

## 3. Non-goals

- **No GPU for CellOracle, SCENIC or CellPhoneDB** — they have no torch backend.
  Making them "GPU-aware" would be theatre.
- **No change to the CPU path's pins.** The seven existing environments and their
  pins are untouched, so no existing user rebuilds anything.
- **No automatic upgrade of the scvi family.** This spec pins what this machine
  already resolves; bumping scvi-tools for its own sake is a separate decision.
- **No multi-GPU / distributed support.** One device, chosen explicitly.
- **No change to `RunAIAnalysis()`'s action surface.** AI actions gain `device`
  in a later milestone, once the underlying functions have it.

## 4. Acceptance criteria

- **A1** `sclet_gpu_status()` reports `available = FALSE` on a machine with no
  NVIDIA driver or no built environment, without erroring and without building
  anything.
- **A2** On this machine: `sclet_gpu_status()` reports the RTX 4090, torch
  `2.13.0+cu130`, and `cuda.is_available() = TRUE`.
- **A3** `RunIntegration(sce, method = "scVI", ..., device = "cpu")` produces
  byte-identical results to the same call without `device` (the default path is
  genuinely unchanged).
- **A4** `device = "cuda"` on a machine without a GPU raises
  `sclet_gpu_unavailable` naming `sclet_gpu_status()`, not a Python error.
- **A5** `device = "auto"` on a CPU-only machine runs on CPU and says so once.
- **A6** `options(sclet.python.device = "cuda")` and
  `SCLET_PYTHON_DEVICE=cuda` both reach the resolver, with the explicit argument
  still winning.
- **A7** The GPU path actually uses the GPU: a callback that records
  `torch.cuda.is_available()`/`torch.cuda.device_count()` inside the worker
  reports what the host reported (pinned by a test that runs only when a GPU is
  present).
- **A8** No existing test changes expectation; the CPU default is not disturbed.
- **A9** A real scVI run on a small object completes on GPU and returns a latent
  representation of the right shape (this machine; gated on A2).

## 5. Tests

New `tests/testthat/test-gpu-backends.R`:

| case | asserts |
|---|---|
| status without env | A1 |
| status on this machine | A2 (`skip_if_not`) |
| cpu default unchanged | A3 |
| cuda without gpu | A4 (simulated by pointing the resolver at an empty `SCLET_ENV_ROOT`) |
| auto fallback | A5 (same simulation) |
| precedence | A6 |
| gpu actually used | A7 (`skip_if_not`) |
| real scVI on GPU | A9 (`skip_if_not`, `skip_if_not_installed("zellkonverter")`) |

The simulation matters: a CI runner has no GPU, and the CPU-fallback behaviour is
what most users will hit, so it must be tested on a machine with no card by
pointing `SCLET_ENV_ROOT` at an empty directory rather than by mocking `torch`.

## 6. Risks

- **3 GB download on first GPU use.** Mitigated by laziness and by documenting it
  in the error/message the user sees the first time.
- **CUDA wheel availability.** `pytorch-cuda` wheels are large and hosted on
  PyPI; a network that cannot fetch them must produce the existing typed
  environment error, not a silent CPU fallback that hides a requested GPU.
- **Driver/CUDA mismatch.** A venv built for cu130 against an older driver fails
  at `import torch`. The typed error must say "driver too old for this CUDA
  build" and name `nvidia-smi`, because that is the actual cause and it is not a
  sclet bug.
- **scvi-tools 1.5 vs 1.4.3 semantics.** The GPU path runs a newer scvi-tools than
  the CPU path, so the two can differ in defaults. A3's byte-identical test is
  CPU-vs-CPU precisely so this is not hidden; the GPU-vs-CPU comparison is
  recorded in NEWS as "not guaranteed identical".
- **Existing pins.** Nothing in `R/basilisk.R`'s seven environments changes, so no
  spec hash moves and no existing environment is rebuilt.

## 7. Rollout

1. `sclet_gpu_env` + `sclet_gpu_status()` + the resolver, with tests A1–A6.
2. `RunIntegration(device = )` GPU path, A7 + A9.
3. `RunSpatialDeconvolution(device = )` GPU path.
4. NEWS + `python-notes.md` + the book's integration chapter.
5. (Later milestone) `device` on the AI actions.

## 8. Open questions for review

1. **Pin the newest scvi family (1.5.1) or match the CPU env (1.4.3)?** Newest is
   what this machine resolves and gets the CUDA build; matching the CPU env keeps
   the two paths semantically identical but may not have a cu130 wheel.
   Recommend: newest, and state the divergence in NEWS.
2. **RegVelo GPU**: its env pins `pytorch=2.5.1` + `scvi-tools=1.1.6`, a different
   generation. A GPU variant means a second regvelo environment (~3 GB). Recommend
   deferring it and saying so, rather than shipping half of it.
3. **`auto` as a default for `RunIntegration` only?** Integration is the one
   method where the speed difference is largest (scVI trains for minutes). Recommend
   no: the default stays CPU everywhere, for the reason in the header.
4. Should `sclet_gpu_status()` be exported (it is a user-facing diagnostic, like
   `sclet_env_path()`) — recommend yes.
