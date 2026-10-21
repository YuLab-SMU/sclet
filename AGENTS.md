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
- **`Collate` must list every file in `R/`.** The committed `DESCRIPTION` listed
  84 of 95 and omitted all eleven `ai-registry-*.R`; the package could not build
  (`R CMD check` dies at install with "files ... missing from 'Collate'"), and it
  had gone unnoticed because `devtools::check()` regenerates docs first and that
  pass re-added them. **`make check` is now `devtools::check(document = FALSE)`**,
  so the tree stays clean; document with `make rd`, which `make all` already runs
  before `check`.
- **An exported object with roxygen but no committed Rd breaks the build.**
  `ValidateAIEvidenceRefs` and `sclet_ai_evidence_independence` were in that
  state. Public Rd is tracked; only internal `man/sclet_ai_*.Rd` is deliberately
  untracked, and that rule does not cover exported functions.

### Basilisk / Python

Read [`.dev/python-notes.md`](.dev/python-notes.md) before touching
`R/basilisk.R`, a Python-backed action, or a Python env. Short version: pins are
paid for by failures (don't tidy a version number), `full.activation = TRUE` on
conda envs, callbacks must be package-scope, `convert = FALSE` when threading
Python objects, and the env path is keyed by the package version so a
DESCRIPTION bump forces a ~10 min rebuild.

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

- **Design semantics may be exposed; analysis results may not.** The privacy
  layer's job is not "no colData values leave" — a `group` contrast cannot be
  planned without the levels of the design column, and withholding them is what
  made one of twelve same-tree trials die with `contrast_invalid` (2026-10-21).
  The line is: **`ActiveIdent` (the clustering) is the result the analysis is
  meant to discover** — it goes out only as `group_N` codes; **a design column
  (`condition`, `Treatment`, `donor`) is the user's own experimental setup** —
  its levels must be present, or the planner has to invent them. Concretely,
  `capabilities$contrast_levels` lists design levels and **must exclude the
  active identity column by name**, and `capabilities$identity_levels` stays
  anonymised. Two traps found while implementing it: `levels` must be a
  character *vector* (an unnamed list of scalars is redacted as "unnamed nested
  values are not bounded"), and `note` is on the privacy deny list, so
  explanatory text belongs in `guidance`.
- **Test fixtures for label-leak tests must not collide with design columns.**
  `test-ai-ident-groups.R` used `treat1`/`treat2` as cluster names, which are
  also `mockSCE()`'s `Treatment` levels — so once design levels were legitimately
  exposed, the "the string appears nowhere" assertion failed on the *design*
  value, not on a leak.

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
- **A previous round is TWO records, never one.** The plan comes from the newest
  `ai_execution` record's `params$steps`; the report from the newest `ai_report`
  record, whose `warnings` live in `artifacts$warnings` (a list of character
  vectors) and whose cited decisions are `summary$cited_user_decisions`. Reading
  either half from the other place yields `not_comparable` on every real object —
  which is indistinguishable from "the judge was never run".
- **`RunDecontX()` writes its record under the fixed id `"decontX"`**, and
  `RecordAIEvidence(source = )` requires a record id that exists. An evidence
  node's name and the record id it derives from are two different things: the
  node can vary by `name`, the source cannot.
- **celda's decontX is deterministic** on a fixed input (identical contamination
  vectors, max abs diff 0), so an identical re-run is an idempotent no-op and that
  is correct — do not pin "a re-run conflicts" for it. A *changed* input under the
  same node name is the case that must conflict.
- **Evidence values admit no free-text strings** — the privacy layer rejects them
  as unbounded. A structural fact such as a written assay name belongs in the
  record's methods section, not in a finding.

- **A pin test that reads files by relative path is testing the directory, not
  the code.** `list.files("../../R")` passes under `devtools::test()` and fails
  under `R CMD check`, where `R/` does not exist beside the tests (13 assertions
  failed this way, 2026-10-21). Assert against the live function instead — a
  prerequisites call that returns the gate name is both location-independent and
  a stronger check than scanning source text. And `devtools::test()` passing is
  not `make check` passing: any change that adds an export or a parameter needs
  the outer loop once, because codoc and Rd checks only run there.
