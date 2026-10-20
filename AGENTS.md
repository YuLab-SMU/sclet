# AGENTS.md — sclet workspace

Single entry point for the `sclet` workspace. **Each directory is its own
checkout of its own branch** — entering a directory means being on that
branch; a file in one directory is not the same version as in another.

| Directory | Branch | Contents |
|---|---|---|
| `.` (top) | `devel` | sclet R package source (`R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `tests/`) |
| `bookdown/` | `docs` | Bookdown source (`*.Rmd`, `_bookdown.yml`, `references.bib`) |
| `bookdown/docs/` | `gh-pages` | Published site (built `.html`; push to publish) |

`bookdown/` is git-ignored at the top level. Package work happens here on
`devel`; book edits in `bookdown/`; never run bookdown builds from the wrong
directory, and never expect a commit in one checkout to appear in another
without rebuilding/pushing.

```bash
git -C .         branch --show-current   # → devel
git -C bookdown  branch --show-current   # → docs
```

## R package development

- **Run checks through the Makefile**: `make check` for routine validation,
  `make check2` for a tarball check, `make clean` for `*.Rcheck/`. Run the
  focused `testthat` filter for the changed module before a broad check.
- **Sparse assays break reshape/plot code.** Never transpose the full assay;
  materialize just the requested features:
  `as.data.frame(t(as.matrix(assay(sce, a)[features, , drop = FALSE])))`.
- **Keep zero-expression cells in plot/fit data** — sparse `summary()`
  triplets omit structural zeros and bias curves upward.
- **Prefer positional alignment over key joins**: `cbind` + melt is ~20x
  faster than `merge(by = "cell")` on large cell counts.

### Test mocks and roxygen2

- **Mock the shared wrapper, don't build envs**: with_mocked_bindings on
  `sclet_basilisk_run`; with `export_all = FALSE` internal `sclet_*` calls
  fail — harness misconfiguration, not a bug.
- **`devtools::document()` rewrites everything** (NAMESPACE reformat/imports
  dropped, roxygen2 version bump, internal Rd regenerated → codoc/link
  failures in `make check`). Working procedure: run document(); restore
  NAMESPACE/DESCRIPTION; delete untracked `man/sclet_ai_*.Rd` (internal Rd
  are deliberately untracked); hand-add only new exports to NAMESPACE and
  hand-add new `Collate:` entries in DESCRIPTION (missing Collate = install
  failure).
- **`devtools::load_all()` and `library(sclet)` are two different trees** —
  verify against `library()` only after `R CMD INSTALL .`.

### Basilisk / Python

- **`packages` are pip-installed; `channels` is a no-op** (basilisk ≥ 1.24).
  Non-PyPI or conda-only packages go in `pip = c(...)`. Introspect
  `env@packages` to assert pins in tests.
- **Errors surface wrapped**: read the inner `error:` of
  `checkForRemoteErrors` and the full traceback via
  `reticulate::py_last_error()`.
- **Use `full.activation = TRUE`** for conda envs (libstdc++ ordering →
  GLIBCXX errors otherwise); keep callbacks package-scope (serialized to the
  worker; closures capturing S4 objects break); `basiliskRun` may be called
  repeatedly on the same proc.
- **Use `convert = FALSE`** when threading Python objects through multi-step
  libraries; only `py_to_r()` the final result.
- **Shims in the shared wrapper** (guarded by `hasattr`): pandas
  `get_values` → `to_numpy`; `np.object` → `object`. Wrap import and env
  setup failures as typed errors (`sclet_glibcxx_not_found`,
  `sclet_env_setup_failed`) with actionable repair text.
- **pySCENIC**: `prune2df(rnkdbs, modules, motif_annotations, ...)`; pass
  `client_or_address = "custom_multiprocessing"` to prune2df; pin the full
  set (`setuptools==80.9.0`, `arboreto==0.1.6`, `dask==2024.2.1`,
  `distributed==2024.2.1`, `dask-expr==0.5.3`, `pyarrow==20.0.0`).

### Python environments on this machine

Find an existing env before concluding a dependency is missing. `uv`
(`~/.local/bin/uv`, no config/proxy overrides needed — use the global proxy
rule `127.0.0.1:7897` on timeouts) creates throwaway venvs: `uv venv
--python 3.12`, then `uv pip install` (no activation; use `./.venv/bin/python`).
Hardlink warnings under `/tmp` are harmless.

sclet's basilisk envs live under
`~/.cache/R/basilisk/<ver>/<pkg>/<pkg-ver>/<envname>/` (declared in
`R/basilisk.R`; seven envs, see pins there). Built lazily — a missing
directory means "not built yet". Path is keyed by package version, so a
DESCRIPTION bump forces a rebuild (long first run; not a network failure).
`obtainEnvironmentPath()` CREATES the env — prefer reading the path directly.
The velocity tests execute in `velociraptor/1.22.0/env/`, not sclet's own.

## AI-layer rules (from the 2026-10-16/17 benchmark sweep; `../sclet-benchmark-ai`)

Each rule below was discovered by a real failure and will bite again if
ignored.

### Standing workflow

- **Records and coverage update in the same commit** — a coverage table that
  under-reports measured work invites re-running done work and hides gaps.
- **A fix without a pinned test/benchmark case is a claim.**
- **When a run contradicts a record, re-grade the record** honestly
  (e.g. the engine's "0 evidence nodes" was right only by accident).

### Running AIAction benchmarks (each rule cost a billed call)

- **`ExecuteAIPlan()` must receive `validation =` explicitly** — passing only
  `confirmation =` makes the engine re-validate and mint a fresh token that
  can never match. After `AIPlanAnalysis()`, take both `validation=` and
  `confirmation=` from ONE fresh strict `ValidateAIPlan()` result; the token
  is never on `plan$validation`.
- **`record_result()` gets fully named flat lists**: nested `c()` loses
  fields (jsonlite invents keys), `do.call(c, ...)` stringifies everything;
  pass each named list as its own argument to `c()`. Unnamed condition lists
  make loops iterate zero times while the runner prints ALL DONE.
- **Pre-flight for billed runs**: `compgen -G ~/Rlib/00LOCK-*` (library
  quiet), `timeout --signal=KILL`, out-of-process heartbeat on record-count
  (a silent no-op is visible there, invisible in the log).

### Measured engine facts

- **ScaleData sparse pathology** (fixed `cbb3590`): row-subset assignment on
  a dgCMatrix hung 20 minutes, and row sds were ONE scalar (the "scaled"
  layer had never been row-scaled). Now rbind2-spliced, true sparse row
  variance, 3.6 s. The model was right to write
  `features = ${step.output.variable_features}` — hang on a benchmark run
  means a data-shape-dependent handler pathology, not the provider.
- **Rare-population seeding** (fixed `04f00dc`): was one label per dense
  cell (2419 populations / 2700 cells). Absurd-on-its-face results mean
  check the algorithm before writing anything up.
- **Message–artifact disagreement is a defect class** (5th instance): assert
  the artifact, not the printed message.
- **Empty plans once validated and ran to `status="completed"`** —
  downstream consumers of `status` should treat zero-action executions as
  suspicious.
- **velociraptor `dimred=` passthrough has an upstream shape bug** (scvelo
  recomputes its own PCA; the prerequisite gate still guards). Upstream
  issue — don't re-debug.
- **Synthetic velocity fixtures need amplitude** (amplitude 8 → zero-edge
  velocity graph, an opaque concat error). Working shape documented in the
  benchmark repo's `R/fixtures_velocity.R`.

### AI-layer architecture

- **The export surface is the external-agent contract** — when adding an
  AI-layer function, ask whether the sclet-driver path (external agent, no
  aisdk) needs it; if yes, export with a documented example.
- **Skills (`inst/skills/`) are injected domain knowledge; the engine stays
  the authority.** Knowledge-type failures (ordering, parameter vocabulary)
  go into a skill, not the inline prompt. Gated by
  `options(sclet.ai.use_skills, default TRUE)`.
- **One plan is one round; starting the next round is the caller's job.**
  `ExecuteAIPlan()` walks `plan$actions` in order with per-step retry —
  nothing inside execution rewrites the later steps, and the repair rounds in
  `AIPlanAnalysis()` fire *before* execution against validation errors and
  goal coverage, never against results. A multi-stage analysis is therefore a
  sequence of plans, each validated against the fingerprint of the object it
  runs on; recording a round changes that fingerprint, so the previous plan
  is refused afterwards by design. The fix for "it doesn't adapt" is a new
  round, not conditional steps inside one plan.
- **Test mocks must use the real schema keys** (`proposed_actions`, not
  `actions`) — a wrong-key mock feeds a normalizer that silently ignores it,
  and the test ends up asserting the bug.
