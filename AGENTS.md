# AGENTS.md — sclet workspace

Each directory is its own checkout of its own branch: entering a directory means
being on that branch, and a file in one directory is not the same version as in
another.

| Directory | Branch | Contents |
|---|---|---|
| `.` (top) | `devel` | R package source (`R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `tests/`) |
| `bookdown/` | `docs` | Bookdown source (`*.Rmd`, `_bookdown.yml`, `references.bib`) |
| `bookdown/docs/` | `gh-pages` | Published site (built `.html`; push to publish) |

`bookdown/` is git-ignored at the top level. Package work happens on `devel`;
book edits in `bookdown/`; never run bookdown builds from the wrong directory.

## R package development

`make check` (full), `make check2` (tarball), `make clean` (`*.Rcheck/`).

**Three tiers — `make check` is the outer loop, not the inner one.**

```bash
# 1) ~15 s, after every edit
Rscript -e 'devtools::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-x.R")'
# 2) ~8 min, when a change is done (filter = "^ai-" is 463 s for one area)
Rscript -e 'devtools::test()'
# 3) ~16 min, milestone close-out / before commit only
make check
```

No single test file dominates: all 81 files sum to ~8 min, the slowest is 57 s,
and `test-malignancy-cnv.R` is 29 s (its cases mock `inferCNV`; nothing runs
JAGS). Do not add `SLOW_TESTS` gates — filter by area if you need it faster.
**A timing claim in a doc is a claim: re-measure it before repeating it.**

- **Long runs go through the job system, not `nohup ... &`** — the shell returns
  in seconds, the harness watches that dead shell, and the completion notice
  never fires. Use `run_in_background: true` then `job_output(wait: true)`.
- **Sparse assays break reshape/plot code.** Never transpose the full assay;
  materialize just the requested features:
  `as.data.frame(t(as.matrix(assay(sce, a)[features, , drop = FALSE])))`.
- **Keep zero-expression cells in plot/fit data** — sparse `summary()` triplets
  omit structural zeros and bias curves upward.
- **Prefer positional alignment over key joins**: `cbind` + melt is ~20x faster
  than `merge(by = "cell")` on large cell counts.

### Test mocks and roxygen2

- **Mock the shared wrapper, don't build envs**: `with_mocked_bindings` on
  `sclet_basilisk_run`; with `export_all = FALSE` internal `sclet_*` calls fail
  — harness misconfiguration, not a bug.
- **`devtools::document()` rewrites everything** (NAMESPACE reformat, dropped
  imports, roxygen2 version bump, regenerated internal Rd → codoc/link
  failures). Working procedure: run document(); restore NAMESPACE/DESCRIPTION;
  delete untracked `man/sclet_ai_*.Rd` (internal Rd are deliberately
  untracked); hand-add new exports to NAMESPACE and new `Collate:` entries to
  DESCRIPTION (missing Collate = install failure).
- **`devtools::load_all()` and `library(sclet)` are two different trees** —
  verify against `library()` only after `R CMD INSTALL .`.

### Basilisk / Python

- **`packages` are pip-installed; `channels` is a no-op** (basilisk ≥ 1.24).
  Non-PyPI or conda-only packages go in `pip = c(...)`. Assert pins by
  introspecting `env@packages`.
- **Errors surface wrapped**: read the inner `error:` of
  `checkForRemoteErrors`; full traceback via `reticulate::py_last_error()`.
- **`full.activation = TRUE`** for conda envs (libstdc++ ordering → GLIBCXX
  errors otherwise); callbacks must be package-scope (serialized to the worker;
  closures capturing S4 objects break); `basiliskRun` may be called repeatedly
  on one proc.
- **`convert = FALSE`** when threading Python objects through multi-step
  libraries; `py_to_r()` only the final result.
- **Shims in the shared wrapper** (guarded by `hasattr`): pandas `get_values` →
  `to_numpy`; `np.object` → `object`. Wrap import/env failures as typed errors
  (`sclet_glibcxx_not_found`, `sclet_env_setup_failed`) with repair text.
- **pySCENIC**: `prune2df(rnkdbs, modules, motif_annotations, ...)` with
  `client_or_address = "custom_multiprocessing"`, and pin the exact set:
  `setuptools==80.9.0`, `arboreto==0.1.6`, `dask==2024.2.1`,
  `distributed==2024.2.1`, `dask-expr==0.5.3`, `pyarrow==20.0.0`.

### Python environments on this machine

Find an existing env before concluding a dependency is missing. `uv` creates
throwaway venvs (`uv venv --python 3.12` then `uv pip install`; use
`./.venv/bin/python`, no activation). Hardlink warnings under `/tmp` are
harmless.

sclet's basilisk envs live under
`~/.cache/R/basilisk/<ver>/<pkg>/<pkg-ver>/<envname>/` (declared in
`R/basilisk.R`; seven envs, pins there). Built lazily — a missing directory
means "not built yet". The path is keyed by package version, so a DESCRIPTION
bump forces a rebuild (long first run, not a network failure).
`obtainEnvironmentPath()` CREATES the env — prefer reading the path directly.
Velocity tests execute in `velociraptor/1.22.0/env/`, not sclet's own.

## AI-layer rules

Every rule below was found by a real failure and will bite again if ignored.

### Standing workflow

- **Records and coverage update in the same commit** — an under-reporting
  coverage table invites re-running done work and hides gaps.
- **A fix without a pinned test/benchmark case is a claim.**
- **When a run contradicts a record, re-grade the record** honestly.

### Running AIAction benchmarks (each rule cost a billed call)

- **`ExecuteAIPlan()` must receive `validation =` explicitly** — passing only
  `confirmation =` makes the engine re-validate and mint a token that can never
  match. After `AIPlanAnalysis()`, take both from ONE fresh strict
  `ValidateAIPlan()` result; the token is never on `plan$validation`.
- **`record_result()` gets fully named flat lists**: nested `c()` loses fields
  (jsonlite invents keys), `do.call(c, ...)` stringifies everything, and
  unnamed condition lists make loops iterate zero times while the runner prints
  ALL DONE.
- **Pre-flight for billed runs**: `compgen -G ~/Rlib/00LOCK-*` (library quiet),
  `timeout --signal=KILL`, and an out-of-process heartbeat on record count — a
  silent no-op is visible there and invisible in the log.

### Measured engine facts

- **Row-subset assignment on a `dgCMatrix` can hang or lie** (hung 20 min, and
  the "scaled" layer held one scalar sd for the whole matrix). A hang on a
  benchmark run means a data-shape-dependent handler pathology, not the provider.
- **Absurd-on-its-face results mean check the algorithm before writing anything
  up** (rare-population seeding once produced one label per cell).
- **Message–artifact disagreement is a defect class**: assert the artifact, not
  the printed message.
- **A write path that swallows its error yields a citable artifact describing an
  older run.** Evidence ids derive from the analysis name, so a re-run collides;
  `sclet_ai_record_evidence_for_run()` is the one recorder — identical values are
  a no-op, different values without `replace = TRUE` raise
  `evidence_id_conflict`. Never wrap a write in `tryCatch(error = function(e)
  object)`.
- **Empty plans can validate and run to `status = "completed"`** — treat
  zero-action executions as suspicious.
- **velociraptor's `dimred=` passthrough has an upstream shape bug** (scvelo
  recomputes its own PCA; the prerequisite gate still guards). Don't re-debug.
- **Synthetic velocity fixtures need amplitude** (amplitude 8 gives a zero-edge
  velocity graph and an opaque concat error). Working shape:
  `../sclet-benchmark-ai/R/fixtures_velocity.R`.

### AI-layer architecture

- **The export surface is the external-agent contract** — when adding an
  AI-layer function, ask whether the sclet-driver path (external agent, no
  aisdk) needs it; if yes, export with a documented example.
- **Skills (`inst/skills/`) are injected domain knowledge; the engine stays the
  authority.** Knowledge-type failures (ordering, parameter vocabulary) go into
  a skill, not the inline prompt. Gated by `options(sclet.ai.use_skills)`.
- **One plan is one round; starting the next round is the caller's job.**
  `ExecuteAIPlan()` walks `plan$actions` in order with per-step retry — nothing
  in execution rewrites later steps, and `AIPlanAnalysis()` repairs fire before
  execution against validation errors and goal coverage, never against results.
  A multi-stage analysis is therefore a sequence of plans, each validated against
  the fingerprint of the object it runs on; recording a round changes that
  fingerprint, so the previous plan is refused afterwards by design. The fix for
  "it doesn't adapt" is a new round, not conditional steps.
- **Test mocks must use the real schema keys** (`proposed_actions`, not
  `actions`) — a wrong-key mock feeds a normalizer that silently ignores it and
  the test ends up asserting the bug.
