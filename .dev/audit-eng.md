# Engineering hygiene audit — sclet (`devel`, 2026-10-22)

Independent, fresh audit. Scope: (1) test suite as a suite, (2) pin-test trap,
(3) reproducibility, (4) packaging from a clean checkout.
Read-only apart from this file. Every finding carries a file:line and a reproduce step.

Status: **COMPLETE** (2026-10-22). All four areas covered; every claim below is
either reproduced by a command shown inline or explicitly marked *not verified*.

---

## 0. Method and what was measured

Read-only. The only file written is this one. Nothing under `R/`, `man/`,
`tests/`, `NAMESPACE` or `DESCRIPTION` was modified, and no build was run inside
the working tree — every build/check was done on copies under `/tmp`.

Tiers actually used (the cheap ones, as instructed):

| What | How | Cost |
|---|---|---|
| Export→test map | `getParseData()` over 93 test files, matched against `NAMESPACE` exports that are also defined in `R/` | seconds |
| Mock census | `grep` for `with_mocked_bindings` / `local_mocked_bindings` per file | seconds |
| Clean-checkout build | `git archive HEAD \| tar -x -C /tmp/sclet-clean`, then `R CMD build` | **2.6 s**, exit 0 |
| Clean-checkout check | `R CMD check --no-tests --no-examples --no-vignettes --no-manual` on that tarball | **4 m 56 s** |
| Reproducibility spot-tests | ad-hoc `Rscript` under `/tmp`, `devtools::load_all()` once per script (~12 s each) | ~40 s/script |

`devtools::test()` (15 min) and `make check` (25 min) were **not** run. One
`R CMD check` on a tarball was run because the packaging question needed a real
answer rather than a reading of the DESCRIPTION.

### A trap worth recording: the stale tarball

`/tmp/sclet_1.1.0.tar.gz` was already present, dated **2026-10-21 21:12**. My
first `R CMD check` silently used it instead of my own build (my build lands in
`/tmp/sclet-clean/`). It reported **3 WARNINGs** — duplicated Rd aliases,
missing Rd cross-references, undocumented arguments. All of them are gone at
`HEAD`: they were fixed by commits `74cf668` (00:59) and `014766a` (01:36) later
that night, and `man/advanced-ai-diagnostics.Rd` (the duplicated-alias source)
no longer exists.

Checking the current tree gives **`Status: 1 NOTE`**. Every Rd finding above was
a false positive about the present code. Reproduce the distinction:

```bash
ls -l /tmp/sclet_1.1.0.tar.gz                    # stale, 728 K, built 21:12
tar tzf /tmp/sclet_1.1.0.tar.gz | grep advanced-ai-diagnostics   # present -> stale
git archive HEAD | tar -x -C /tmp/sclet-clean && (cd /tmp/sclet-clean && R CMD build .)
tar tzf /tmp/sclet-clean/sclet_1.1.0.tar.gz | grep advanced-ai-diagnostics  # absent -> current
```

**Lesson for the next auditor: check the mtime and contents of any tarball you
find in `/tmp` before you believe a check log.** A 13-hour-old tarball reports a
package that no longer exists.

---

## 1. The test suite as a suite

Corpus measured: **93 `test-*.R` files** (the brief says 91 — two were added since),
**849 `test_that()` blocks**, 95 files under `tests/`. No `helper*.R` file exists,
so every file builds its own fixtures.

### 1.1 Exported functions with no test

Method: take the `export(...)` list from `NAMESPACE` (281 names), drop the 13 that
are re-exports not defined in `R/`, leaving **268 exported functions of sclet's own
code**. Then parse all 93 test files with `getParseData()` and collect every
`SYMBOL_FUNCTION_CALL` token. This counts both `RunPCA(...)` and `sclet::RunPCA(...)`
and — unlike `grep` — does **not** count a name that appears only in a comment or a
string. Reproduce with the script in §0's table.

**Result: 105 of 268 exported functions (39%) are never called by name in any
test.** Split by kind:

| kind | untested |
|---|---|
| `Run*` analysis modules | **28** |
| plotting (`plot_*`, `DimPlot`, `VlnPlot`, …) | **40** |
| accessors (`get_*` / `has_*`) | 19 |
| other | 18 |

The accessors and plots are individually cheap; the 28 modules are not. The full
list is reproducible via the script; the ones that matter most:

**Most surprising gaps — functions that compute a statistic:**

| function | why it matters | defined at |
|---|---|---|
| `Read10X` | the first call in the README; parses 10x output | `R/io.R:18` |
| `subset_feature` | first filtering step in the README | `R/subset.R:26` |
| `PercentageFeatureSet` | computes per-cell QC metric | `R/preprocessing.R:36` |
| `QCMetrics` | computes the QC table everything downstream uses | `R/preprocessing.R:14` |
| `FindMarkers_pseudobulk` | pseudobulk aggregation + test | `R/pseudobulk.R:85` |
| `refit_milo` | differential abundance refit | `R/refit_milo.R:41` |
| `RunPhenotypeAssociation` | association test | `R/phenotype.R:23` |
| `RunGeneSetScoring` | module scoring | `R/geneset_scoring.R` |
| `RunGRN`, `RunRegVelo` | regulon / regulatory velocity | `R/scenic.R:419`, `R/regvelo.R:42` |
| `RunMilo`, `RunVelocity`, `RunCellRank`, `RunCellFate` | mainline statistics | `R/runMilo.R:72`, `R/velocity.R:25`, `R/cellrank.R:30` |
| `AuditAnalysisChain`, `PipelineSummary` | summarise an analysis for a reader | `R/copilot.R:77`, `R/pipeline_summary.R:12` |
| `WriteH5AD`, `as.SCE` | round-trip; silently losing state is the risk | `R/interop.R:40`, `R/interop.R:95` |

The first four are the striking ones: they are exactly the code the README tells a
newcomer to run in step 1, and nothing pins them. I ran them by hand (§4.3) and
they work — but "works today on my machine" is what the suite exists to replace.

**Honest caveat on the count.** "Never called by name" is not "never executed". A
few of the 28 modules are reached through the AI action registry with a **mocked
handler**, which tests the wiring, not the computation:

```
run_milo        -> 3 test files        run_velocity   -> 7
run_decontx     -> 2                   run_cellrank   -> 0
run_spatial_niche 0   run_imputation 0   run_geneset_scoring 0
run_phenotype  0       run_regvelo    0   run_slingshot       0
```

So the true figure is "≈100 functions with no test that executes them", not
"105 functions with no test at all". Both numbers say the same thing.

### 1.2 Mocking versus real logic, per file

The suite is unusual: **it barely mocks.** Only **12 of 93 files** call
`with_mocked_bindings()` / `local_mocked_bindings()`, **25 calls in total**:

| mocks | file |
|---|---|
| 5 | `test-ai-fate.R` |
| 4 | `test-ai-execution.R` |
| 3 | `test-ai-communication.R` |
| 2 | `test-knn-compat.R`, `test-ai-skills.R`, `test-ai-interpretation-gaps.R`, `test-ai-evidence-rerun.R` |
| 1 | `test-scenic.R`, `test-ai-single-cell-chain.R`, `test-ai-routes.R`, `test-ai-plan-error-guidance.R`, `test-ai-evidence-linked-report.R` |

The other 81 files run real code on `scuttle::mockSCE()` fixtures. That is the
right choice for fidelity, and it is why the suite takes 15 minutes.

**But the strategy that replaces mocking is `skip_if_not_installed()`, and that is
the suite's main structural weakness.** 211 guards across 35 files:

```
scuttle 40   aisdk 22   velociraptor 19   slingshot 16   BiocNeighbors 14
scran 13     scater 12  igraph 11         clusterProfiler 10
SingleCellExperiment 8   org.Hs.eg.db 8   ggplot2 7   basilisk 6   batchelor 5
```

Consequence: **how much of the suite actually runs depends on what is installed.**
The CI matrix runs with `_R_CHECK_FORCE_SUGGESTS_ = false`, so a cell that is
missing a Suggests does not fail — it **skips**, and the run is green. `test-state-refactor.R`
alone has 41 skip guards; `test-ai-execution-summary-contract.R` has 19. A file
with 41 guards has 41 ways to report nothing while looking like it passed. That is
not fixable by mocking everything, but the suite should not be described as "93
files of coverage" without that qualifier, and a CI cell that installs the full
Suggests set is what makes the number mean anything.

Note also: `scuttle` (the fixture source) is guarded 40 times, so on a machine
without it nearly the whole suite skips at once — a single missing package can
silently empty the run.

### 1.3 Tests that break under a harmless refactor

Ranked by how much maintenance each costs.

**Class A — breaks on a renamed internal or reworded source (the worst; these
should be rewritten).** Two files, detailed in §2.2:

- `test-ai-evidence-rerun.R:70-72` — 8 assertions grepping `deparse(body(fn))` for
  the string `sclet_ai_record_evidence_for_run`. Renaming an internal breaks 8
  assertions with no behaviour change.
- `test-env-path-version-free.R:91-92` — matches `LD_LIBRARY_PATH=` in the
  deparsed body of `sclet_env_run`. Reformatting the same string breaks it.

**Class B — breaks on a reworded message.** 97 `expect_error()` calls match a
**regex against the message text** (`expect_warning`: 2; `expect_message`: 0, so
the exposure is contained). Most patterns are short and stable; these are full
sentences and will break on any copy-edit:

```
"refusing to record an AI result citing unresolvable evidence"
"refusing to record an over-claimed AI result"
"could not be resolved from any available source"
"geneset_scoring: No gene set scoring results found"
"No Bioconductor-native HVG statistics found"
"proposed_claim_level must be one of"
```

Counter-example of good practice in the same corpus: `"over-claimed|unresolvable|evidence ref"`
is an alternation that survives rewording. `fixed = TRUE` appears 293 times, and
substring matching (as opposed to regex) is already the dominant style — the
long-sentence `expect_error` patterns are the minority that will bite.

**Class C — breaks on a reordered list.** Only 4 assertions across the corpus match
`expect_identical(names(` / `expect_identical(rownames(` and similar. That is a
very low number for a 849-block suite, and it means list ordering is barely pinned.
Good for maintenance, mild loss of contract checking.

**Class D — breaks because it depends on the machine, not the code.** The `/tmp`
probe in `test-ai-round-progression.R:239` (§2.3), plus the 211 skip guards above.

**What does *not* break:** there is no `expect_snapshot()` anywhere (0 uses; the
`tests/testthat/_snaps/` directory exists but is empty, and `.Rbuildignore` excludes
it — harmless). Snapshot tests are the usual source of refactor churn and the
project has avoided them entirely.

### 1.4 Are there test files that are effectively duplicates?

**No. Measured, not guessed.**

- **Identical `test_that()` descriptions in more than one block: 0** out of 849
  blocks. (Parsed with a multi-line-aware regex; a naive line-based regex produces
  a false hit on `get_program(source='` in `test-program-resolution.R`.)
- **Call-set similarity** between every pair of the 93 files (Jaccard over the set
  of function names each file calls). The highest pair is **0.594** —
  `test-ai-plan-stop-conditions.R` vs `test-ai-velocity.R`, sharing 38 names — and
  that overlap is scaffolding (`new_sclet_ai_plan`, `sclet_ai_evidence_get_all`,
  …), not assertions. The next pairs are 0.592 and 0.538. A true duplicate would
  score far higher; nothing here is a copy with a new name.

The file names cluster (`test-ai-plan-*`, `test-ai-velocity-*`) but the contents
do not overlap in what they pin. This is the one area of the suite with no finding.

---

## 2. The pin-test trap

**Verdict: the recorded defect class is gone. Two instances of the same failure
mode survive in a different shape, and they are the maintenance cost.**

### 2.1 The recorded pattern: fixed

AGENTS.md (lines 248-255) records 13 assertions that read `list.files("../../R")`
and passed under `devtools::test()` but failed under `R CMD check`. Grepping the
exact patterns the brief names finds **one hit, and it is a comment describing the
fix**:

```
tests/testthat/test-ai-data-support.R:59:
    # An earlier version scanned list.files("../../R"), which resolves in the
```

```bash
cd tests/testthat
grep -rn 'list\.files("\.\.\|readLines("\.\.\|file\.path("\.\.\|dir("\.\.\|system\.file(' .
# -> exactly 1 line, a comment
```

Independently confirmed by the clean-checkout run: `R CMD check` on a tarball
built from `git archive HEAD` reports

```
* checking dependencies in `tests` ... OK
```

and would have failed those assertions had any survived. No test reads `R/`,
`man/`, or the package source directory by path.

The only textual mentions of `R/*.R` in the corpus are documentation comments
pointing a reader at the implementation (`test-ai-fate.R:52`,
`test-ai-fate-drivers.R:25`, `test-ai-spatial-multimodal.R:31-40`,
`test-ai-communication.R:47`, `test-phenotype-evidence.R:13`,
`test-ai-hvg-reference.R:111`, `test-ai-spatial-multimodal.R:289`,
`test-regvelo.R:259`). Comments are not assertions — harmless.

### 2.2 Same failure mode, new shape (1): asserting on `deparse(body(...))`

Two tests assert on a function's **source text** rather than its behaviour. These
are location-independent, so they survive `R CMD check` — but they break on a
harmless refactor, which is the property that makes the pin-test trap expensive.

**Finding 2.2-a — `tests/testthat/test-env-path-version-free.R:91-92`**

```r
body_text <- paste(deparse(body(sclet_env_run)), collapse = "\n")
expect_match(body_text, "LD_LIBRARY_PATH=", fixed = TRUE)
```

This asserts the *string* `LD_LIBRARY_PATH=` appears in the body of
`sclet_env_run()`. It is a proxy for a real behaviour — the child python must not
inherit R's loader path — but it tests the spelling, not the behaviour. Rewriting
the environment as `c("LD_LIBRARY_PATH" = "")` (identical behaviour), or moving
the env construction into a helper called by `sclet_env_run()`, fails the
assertion while the shipped behaviour is unchanged. The file's own comment
(lines 85-89) says the intent is "Asserted on the constructed argument rather
than by running python" — but it never asserts on the argument; it asserts on the
text that builds the argument. The honest version is to call the small helper and
inspect the returned vector.

**Finding 2.2-b — `tests/testthat/test-ai-evidence-rerun.R:59-74`**

```r
for (name in recorders) {
    fn <- getFromNamespace(name, "sclet")
    body_text <- paste(deparse(body(fn)), collapse = " ")
    expect_true(grepl("sclet_ai_record_evidence_for_run", body_text, fixed = TRUE), ...)
}
```

Eight assertions (one per recorder) that each recorder's body **mentions the name**
`sclet_ai_record_evidence_for_run`. The property the test intends is behavioural —
"every recorder routes through the re-run-aware writer, so none silently drops a
node" — but it is checked by grepping source text. Consequences:

- rename the internal (a harmless refactor) → 8 failures, no behaviour change;
- a recorder that reaches the writer indirectly (`do.call`, a wrapper, a local
  alias) fails the test while behaving correctly;
- a recorder whose body mentions the name in a comment or a dead branch passes
  while doing the wrong thing.

The behavioural test is available and cheap, because the same file already builds
the fixture at line 48: record twice with *different* values and expect
`evidence_id_conflict`, which is the observable consequence of routing through the
shared recorder. The assertion at lines 41-44 already does exactly this for one
recorder. Extending that loop over the eight names would test the property instead
of the spelling.

### 2.3 Same failure mode, new shape (2): a test keyed on an absolute `/tmp` path

**Finding 2.3 — `tests/testthat/test-ai-round-progression.R:239-241`**

```r
skip_if_not(file.exists("/tmp/sclet-agentrun/m7-report.rds"),
    "the M7 object is not present (real-data probe, run locally)")
sce <- readRDS("/tmp/sclet-agentrun/m7-report.rds")
```

This is the pin-test trap inverted: instead of *failing* off the author's machine,
it **silently skips**. It reads a file that is not in the package and cannot be in
a clean checkout, so it asserts nothing in CI, in `R CMD check`, or on any other
machine — while still counting toward the 93 test files. On the one machine where
it runs, it runs against whatever is in `/tmp`, so its verdict is not reproducible
either.

To the author's credit this is deliberate and documented (lines 235-238: "a test
that fails on a fresh clone because /tmp is empty teaches the reader to ignore
failures. The shape it verifies is pinned above by fixtures."). Graded as such: it
is not a defect, it is **coverage that exists only on one laptop**. It should not
be counted in any coverage claim, and the shape it covers — the previous-round
progression over a real M7 object — has no test that runs anywhere else.

### 2.4 The wider silent-skip population

The same "asserts nothing, reports green" mechanism is widespread and unavoidable
(the suite mocks Python backends a lot):

```
skip_if_not_installed  211 occurrences across 35 files
skip_if_not              9
bare skip()              2
```

`test-advanced.R:42` is the one that matters, because the skipped test is the
*only* thing standing between `RunSlingshot` and zero coverage:

```r
test_that("RunSlingshot works", {
    if (requireNamespace("slingshot", quietly = TRUE)) {
        ...
        # sce <- RunSlingshot(sce, cluster.labels=colLabels(sce))
        # expect_true(!is.null(sce@metadata$slingshot))
    }
    # The body above is a hand-run snippet with its assertions commented out;
    # keep it as documentation, but do not report a green test that cannot fail.
    skip("RunSlingshot on a colLabels fixture: hand-verified, not automated yet")
})
```

The comment at lines 40-42 is exactly the right instinct and I am not
second-guessing it. But the effect is that `RunSlingshot` — a load-bearing
trajectory entry point — has **no automated test at all**, and the file's name
("RunSlingshot works") makes the suite look better than it is to anyone reading a
test list rather than the file.

Note also that the skipped body is itself non-reproducible: line 32
`matrix(rpois(1000, 5), ...)` and line 33 `matrix(rnorm(100), ...)` both draw
without a seed. If the skip is ever lifted, the fixture must be seeded first (see
§3).

### 2.5 Tests that assert on file names or paths (checked, clean)

- `test-ai-report-layer.R:157-158` asserts `file.exists(path)` then compares
  `readLines(path)` against what the writer returned — but `path` comes from
  `tempfile()` in the test, so it tests the writer, not the directory. Fine.
- No test asserts a file *name* matches a source file name.
- No test reads `DESCRIPTION`, `NAMESPACE`, `Collate`, or `inst/` by relative path.

---

## 3. Reproducibility

**Verdict: strong where it was recently fixed, with one real gap and one
correct-by-design non-determinism.**

Commit `8db446f` (2026-10-22 00:30) is called "feat(seed): record the seed for
every stochastic path". I tested the claim rather than reading it.

### 3.1 What is seeded, and does it hold?

Nine exported functions take an RNG `seed` formal:

```
FindClusters (clustering.R)          RunUMAP (dimred.R)
RunIntegration (integration.R)       ClassifyMalignantEpithelium (malignancy.R)  <- see 3.3
refit_milo (refit_milo.R)            RunRegVelo (regvelo.R)
RunMilo (runMilo.R)                  RunSlingshot (runSlingshot.R)
RunSCENIC (scenic.R)
```

Four of them actually call `set.seed(seed)` in their body — `clustering.R:127`,
`dimred.R` (via the uwot seed argument), `refit_milo.R:63`, `runMilo.R:108`,
`runSlingshot.R:58`. There is also a seed-safe helper in `state.R:862-890` that
restores the caller's RNG state afterwards, which is the right thing to do and is
the kind of bug (`set.seed()` silently moving the user's stream) that usually goes
unnoticed.

**Tested: two runs with the same `seed` under different global RNG states agree.**

```r
set.seed(2); u1 <- RunUMAP(a, seed = 42)
set.seed(3); u2 <- RunUMAP(a, seed = 42)
identical(reducedDim(u1, "UMAP"), reducedDim(u2, "UMAP"))   # TRUE
```

**Tested: the seed is recorded in the analysis state**, not just used.

```r
recs <- sclet:::sclet_get_state(cl)$states$records
# records: preprocess, reduction, graph, clustering
recs$reduction$umap$params$seed              # 42
recs$clustering$louvain_clusters$params$seed # 123
```

So a downstream reader can recover which seed produced the UMAP and the
clustering. That is the property the commit claims and it is true.

### 3.2 What looks random but is not (a false positive I nearly reported)

`RunPCA()` takes **no** `seed` formal, and two runs under different RNG states
produce `reducedDim()` matrices that are **not** `identical()` — maximum absolute
difference **14.5**. Reporting that as "RunPCA is not reproducible" would have
been wrong. The components differ by **sign only**:

```r
cor(ea[,1], eb[,1])  # -1
cor(ea[,2], eb[,2])  # -1
# align signs of each component, then compare:
max(abs(sweep(ea, 2, sign(diag(cor(ea, eb))), "*") - eb))   # 2.14e-06
```

A residual of 2.1e-06 is floating-point noise; the embeddings are the same up to
the sign ambiguity that every PCA has. **Not a defect.** Worth recording anyway:
a user comparing two PCA runs by `identical()` will see `FALSE`, so if the
package ever documents PCA stability it should say "up to sign".

### 3.3 Finding 3.3 — a genuinely unseeded stochastic path (`malignancy`)

This is the one real reproducibility gap.

`ClassifyMalignantEpithelium()` has a formal called `seed` — but it is **not the
RNG seed**. It is the gene signature produced by `FindMalignantSeedGenes()`
(`R/malignancy.R:312`, `@param seed Either the output of FindMalignantSeedGenes()
or a ...`). There is therefore **no way to seed the randomness** in this function.

The randomness is real and it changes results. `ClassifyMalignantEpithelium()`
calls `sclet_add_module_score()` (defined `R/malignancy.R:828-880`), whose formals
are `expr, genes, bins, n_ctrl, cells, exclude` — **no seed** — and which draws
random control genes:

```r
# R/malignancy.R:866
ctrl <- if (length(pool) > n_ctrl) sample(pool, n_ctrl) else pool
```

`set.seed()` appears **nowhere** in `R/malignancy.R`, and the k-means inside runs
with `nstart` (line 703), which is itself random.

Reproduced on a fixed input:

```r
set.seed(1); s1 <- sclet:::sclet_add_module_score(expr, genes, bins = bins, n_ctrl = 3, cells = colnames(expr))
set.seed(2); s2 <- sclet:::sclet_add_module_score(expr, genes, bins = bins, n_ctrl = 3, cells = colnames(expr))
identical(s1, s2)                 # FALSE
format(max(abs(s1 - s2)), digits = 4)   # 0.3377
```

The same 40x30 matrix, called twice, gives module scores differing by 0.34. The
module score is the input to the malignant/normal split, so the *labels* can move
too — and `ClassifyMalignantEpithelium` is the entry point of the malignancy
mainline.

Why it was missed: the scan for "functions that take a seed" classifies
`ClassifyMalignantEpithelium` as seeded, because the argument is literally named
`seed`. **A name-based seed census is wrong here**; the argument is a signature.
This is worth writing down in AGENTS.md next to the other seed rules.

Fix shape: give `sclet_add_module_score()` and `ClassifyMalignantEpithelium()` a
`random_seed =` (or `seed_rng =`) argument, `set.seed()` around the loop and the
k-means, and record it in `params` like the other paths do. Renaming is safer than
reusing `seed`, since `seed` already means the signature in this API.

### 3.4 Non-determinism that is correct and must not be "fixed"

Two unseeded `sample()` calls in `R/ai-planning.R` are deliberate:

- `R/ai-planning.R:127` — `sclet_ai_plan_id()` builds a unique plan id
  (`plan_<task>_<timestamp>_<8 random chars>`). An id must be unique per plan;
  making it deterministic would collide two plans created in the same second.
- `R/ai-planning.R:1114` — `confirmation_token` is a 32-character random token.
  Making it predictable would defeat the confirmation gate's purpose.

Neither is an analysis result, and neither should be seeded. Recording them here so
a future seed audit does not flag them as defects.

### 3.5 Fixtures that would become non-reproducible if unskipped

`tests/testthat/test-advanced.R:32-33` (inside the `skip()` for `RunSlingshot`)
constructs its fixture with unseeded `rpois()` and `rnorm()`. Harmless today
because the block never runs; it must be seeded before that skip is ever lifted.
See §2.4.

### 3.6 Not verified

The Python-backed stochastic paths (`RunCellRank`, `RunCellFate`, `RunVelocity`,
`RunRegVelo`, `RunDecontX`, `RunSuperCell`) were **not** tested end to end — they
need basilisk environments, which are minutes each and outside the cheap tier. I
read that `RunRegVelo` takes a seed and that the others go through
`sclet_basilisk_run`; I did not verify that a Python-level seed is propagated into
`scvelo`/`cellrank`. Claims about those remain unvalidated.

---

## 4. Packaging as someone else would receive it

### 4.1 DESCRIPTION: used-but-undeclared, declared-but-unused

**Verdict: honest.** This is the one part of the package I could not fault.

Evidence, from the real `R CMD check` on the current tree:

```
* checking package dependencies ... INFO
Packages suggested but not available for checking: 'nichenetr', 'scPAS'
* checking dependencies in R code ... OK
* checking for unstated dependencies in examples ... OK
* checking for unstated dependencies in `tests` ... OK
```

`checking dependencies in R code ... OK` means R CMD check found **no** `::` /
`:::` call to a package absent from DESCRIPTION, and no declared Imports entry
that is never imported. `nichenetr` and `scPAS` are merely *not installed here* —
an INFO, not a defect; both are correctly declared in Suggests, and `scPAS` is
also declared in Enhances.

I checked the reverse direction by hand, because R CMD check does not. Every
declared package is referenced somewhere in `R/`, `tests/` or `inst/`:

```r
# 0 of 56 Suggests and 0 of 27 Imports are never mentioned
```

So there is **no declared-but-unused package**, and no used-but-undeclared
package. Specific things I expected to be wrong and were not:

- `tinytable` (Imports) never appears as `tinytable::`, which looks unused. It is
  used through `importFrom(tinytable, tt, group_tt, style_tt, save_tt)` at
  `R/gene-summary-table.R:14-17` and called by bare name. Correctly declared.
- `parallel::detectCores()` at `R/priority.R:20` is undeclared — but `parallel` is
  a **base** package, so R CMD check does not require declaring it. Same for
  `tools::` (`R/basilisk.R:200`, `R/ai-planning.R:167`).
- `BiocManager::install(...)`, `remotes::install_github(...)`,
  `devtools::test()`, `irlba::`, `uwot::`, `matrixStats::`, `dir.expiry::` all
  appear only inside **error-message strings or comments** — not calls. A naive
  `grep 'pkg::'` reports these as undeclared dependencies (13 false positives);
  parsing the AST for `SYMBOL_PACKAGE` tokens removes them all. Worth knowing
  before anyone "fixes" the DESCRIPTION.

One latent coupling, minor and worth a line rather than a section:
`BiocParallel` is in Suggests but is required at runtime by
`sclet_model_gene_var()`:

```r
# R/preprocessing.R:628-630
if (is.null(BPPARAM)) {
    BPPARAM <- getFromNamespace("SerialParam", "BiocParallel")()
}
```

There is no `requireNamespace("BiocParallel")` guard, and `getFromNamespace()` is
invisible to R CMD check's dependency scan — which is why the check says OK. It
works today only because `scran` (Imports) pulls `BiocParallel` in, so in practice
it is always installed. No action needed, but the declaration understates the
real requirement.

### 4.2 Would it build correctly from a clean checkout?

**Yes. Verified end to end. This is the strongest result in the audit.**

The defect class AGENTS.md records (a `Collate` missing eleven `ai-registry-*.R`
files, and exported objects with no committed Rd) is genuinely fixed: **every
build input is tracked.**

```bash
# what a clean clone would contain
git archive HEAD | tar -x -C /tmp/sclet-clean
ls /tmp/sclet-clean/R/*.R    | wc -l   # 95   == git ls-files R/   | wc -l  (95)
ls /tmp/sclet-clean/man/*.Rd | wc -l   # 279  == git ls-files man/ | wc -l  (279)

# disk vs tracked, per directory: nothing untracked
comm -23 <(ls R/ | sort)   <(git ls-files R/   | sed 's|^R/||'   | sort)   # empty
comm -23 <(ls man/ | sort) <(git ls-files man/ | sed 's|^man/||' | sort)   # empty

# and the working tree is content-identical to HEAD except this audit's own files
diff -rq /tmp/sclet-clean . | grep -v '^Only in \.:'
#   Only in ./.dev: audit-*.md, audit-*.R      (this audit; .Rbuildignore'd via ^\.dev($|/))
#   Only in ./tests/testthat: Rplots.pdf       (stray artifact, ignored by git)
#   Only in ./tests/testthat: _snaps           (empty dir; no snapshot tests exist)
```

`git status --short` is clean apart from the four untracked `.dev/` files this
audit and two sibling audits created. **No untracked file under `R/` or `man/`.**
The `man/sclet_ai_*.Rd` files AGENTS.md describes as "deliberately untracked" no
longer exist as a class — `man/` is 279 files, all tracked, all present in the
tarball.

Build and check from that clean export:

```
R CMD build .                     #  2.6 s, exit 0
R CMD check --no-tests --no-examples --no-vignettes --no-manual  <tarball>   # 4 m 56 s
Status: 1 NOTE
```

The single NOTE is:

```
* checking R code for possible problems ... NOTE
Found the following possibly unsafe calls:
File 'sclet/R/zzz.R':
  unlockBinding("obtainEnvironmentPath", ns)
```

That is the basilisk lookup shim documented in AGENTS.md (lines 140-141), working
as designed. It is a plain `R CMD check` finding, not a BiocCheck one, so it is in
scope — but it is the *deliberate* price of the shim, and I am not proposing to
remove it. Everything else is OK: namespace, S3 methods, replacement functions,
codoc, Rd metadata, Rd cross-references, Rd `\usage`, missing documentation
entries, installability, load/unload.

#### Finding 4.2-a — a stray `Rplots.pdf` is shipped in the tarball (minor, real)

`Rplots.pdf` sits at the package root (4.6 KB, dated from a plotting call). It is
listed in `.gitignore` but **not** in `.Rbuildignore`, so any build run in the
working tree — which is exactly what `make build` and `make check2` do — picks it
up:

```bash
cp Rplots.pdf /tmp/sclet-clean/ && cd /tmp/sclet-clean && R CMD build .
tar tzf sclet_1.1.0.tar.gz | grep Rplots
#   sclet/Rplots.pdf          <-- ships to every downstream user
```

`tests/testthat/Rplots.pdf` is the same artifact one level down and is already
excluded by the `tests/testthat/_snaps`-adjacent rules? No — it is excluded only
because `R CMD build` drops `tests/testthat` extras? It is not excluded either;
it simply is not present in a clean export.

**Fix:** add `^Rplots\.pdf$` to `.Rbuildignore` (one line), or delete the file. A
clean checkout is unaffected, which is why `R CMD check` on the `git archive`
tarball reports `checking top-level files ... OK` — the file is only in the
working tree. This is the one way in which "builds from a clean checkout" and
"builds on the maintainer's machine" differ today.

#### Finding 4.2-b — an in-tree tarball can carry `.git` (cosmetic; mechanism not verified)

A tarball built from the *working tree* rather than a git export can contain a
`sclet/.git` entry, which `R CMD check` then reports as:

```
* checking for hidden files and directories ... NOTE
Found the following hidden files and directories:
  .git
These were most likely included in error.
```

Observed in the stale `/tmp/sclet_1.1.0.tar.gz` (built 2026-10-21 21:12):

```bash
tar tzf /tmp/sclet_1.1.0.tar.gz | grep '^sclet/[.]git'
#   sclet/.git
```

I could **not** reproduce it: a clean export has no `.git`, and both
`R CMD build pkg` and `R CMD build --keep-empty-dirs pkg` on a replica package
excluded `.git` entirely (contents *and* the directory entry). So I record it as
"seen once in an old tarball, not reproducible from a clean checkout" rather than
as a live defect. **Mechanism: not verified.** It matters only because it means the
release artifact a maintainer builds in place is not byte-identical to the one
built from a clean clone — a CI release job built from a fresh clone would not see
it.

**The headline for this section:** a fresh clone builds and checks. Nothing needs
to be recovered from an untracked file. That is a real improvement over the state
AGENTS.md describes, and the `Collate`/Rd defects it records are fixed.

### 4.3 README accuracy

**Verdict: accurate, and it works. A newcomer can get to a running analysis in one
step.**

Verified mechanically:

- **Every function the README names exists and is exported.** Extracting all 67
  `Name()` tokens from `README.md` and intersecting with the 281 `NAMESPACE`
  exports leaves **zero** missing. (One apparent miss, `Workflow()` at line 16, is
  prose — it is the `Run<Mainline>Workflow()` template.)
- **Every documented argument is a real formal.** For `RunBasicWorkflow`,
  `subset_feature`, `DimPlot`, `RunDEtest`, `RunEnrichment`, `EnrichmentPlot`,
  `PercentageFeatureSet`, `Read10X`, `get_detest`, `Status`: all named arguments
  are in `formals()`.
- **The "eleven more" claim is exactly right.** README line 74 names
  `check_de_readiness()`, `check_enrichment_readiness()`, `check_fate_readiness()`
  "and eleven more" = 14. `NAMESPACE` exports **14** `check_*_readiness`
  functions. Someone counted.
- **`README.md` is in sync with `README.Rmd`.** Rendering is 12 s after the source
  (2026-10-21 10:45:18 → 10:45:31) and the only diffs are pandoc reflow
  (` ```r` → ` ``` r`, line wrapping). No unrendered ` ```{r` chunks leaked into
  the `.md`. Both are `.Rbuildignore`d, so neither ships — correct.

Verified by running it (script at `/tmp/readme_quickstart.R`): because there is no
10x directory on this machine, step 1 is substituted with `mockSCE()` and the rest
is the README verbatim.

```r
sce <- scuttle::mockSCE(ncells = 200, ngenes = 400)   # stands in for Read10X()
sce[["percent.mt"]] <- PercentageFeatureSet(sce, "^MT-")
sce <- subset_feature(sce, mincell = 3, peek = FALSE)
sce <- RunBasicWorkflow(sce, n_features = 500, n_pcs = 15, cluster_resolution = 0.8)
Status(sce)
p <- DimPlot(sce, reduction = "UMAP")
sce <- RunDEtest(sce, ident.1 = "1", ident.2 = "2")
head(get_detest(sce)$artifacts$result)
```

Output: `ALL QUICKSTART STEPS 1-4 PASSED`, with a real DE table
(`Gene_0037, pval 6.4e-07, padj 1.9e-04`) and `DimPlot` returning a `ggsc` object.
Both `RunBasicWorkflow` and `subset_feature` accepted the README's argument names
on the first try.

One documentation-grade observation, not a defect: `DimPlot()` returns a **`ggsc`**
object, not a plain `ggplot`. A newcomer who then calls
`p + ggplot2::theme_minimal()` gets an error they will not understand from the
README. Worth one line in the README ("plots return `ggsc` objects") rather than a
code change.

**Not verified:** step 5 (`RunEnrichment()` / `EnrichmentPlot()`) was not executed —
it needs `org.Hs.eg.db` and a reachable annotation source, and I did not want a
network call in a read-only audit. Its arguments are real formals (§above), so the
call will at least dispatch.

---

## 5. Scorecard and priorities

| Area | Grade | One line |
|---|---|---|
| 1. Test suite as a suite | **weak** | 849 blocks and no duplicates, but **105/268 (39%) of own exports are never executed by a test**, and 211 `skip_if_not_installed` guards mean a green run can assert very little. |
| 2. Pin-test trap | **strong** | The recorded defect class is gone — verified by `R CMD check` on a clean checkout, not by grep. Two source-text assertions remain and are the main refactor cost. |
| 3. Reproducibility | **adequate** | Seeded paths are genuinely reproducible *and* record their seed in the analysis state; one exported mainline (`malignancy`) has real unseeded randomness. |
| 4. Packaging | **strong** | A clean clone builds in 2.6 s and checks to **1 NOTE**; DESCRIPTION has no undeclared/unused package; the README is accurate and its quickstart runs. Two minor defects (`Rplots.pdf` shipping; an unreproducible `.git` in an old tarball). |

### Top three findings

**1. Two fifths of the exported surface is never executed by a test — including
the README's own first four calls.** `NAMESPACE` (268 own exports) vs the parsed
test corpus: 105 functions are never called by name, and for ≈100 of them no test
executes the computation either (a few have registry-wiring tests with mocked
handlers). The ones that matter:
`Read10X` (`R/io.R:18`), `subset_feature` (`R/subset.R:26`),
`PercentageFeatureSet` (`R/preprocessing.R:36`) and `QCMetrics`
(`R/preprocessing.R:14`) — these are literally
README quickstart steps 1-2, and nothing pins them. Then 28 `Run*` analysis
modules including `RunMilo` (`R/runMilo.R:72`), `RunVelocity` (`R/velocity.R:25`),
`RunCellRank` (`R/cellrank.R:30`), the pseudobulk path
`FindMarkers_pseudobulk` (`R/pseudobulk.R:85`), `refit_milo`
(`R/refit_milo.R:41`) and `RunPhenotypeAssociation` (`R/phenotype.R:23`). Reproduce with
the AST script in §0. This is the highest-value work: one test per `Run*` module
that reaches its readiness gate and its record contract would move the number more
than any other single change.

**2. `ClassifyMalignantEpithelium()` is stochastic with no way to seed it, because
its `seed` argument is a gene signature.** `R/malignancy.R:866`
(`ctrl <- if (length(pool) > n_ctrl) sample(pool, n_ctrl) else pool`) lives in
`sclet_add_module_score()` (`R/malignancy.R:828-880`), whose formals are
`expr, genes, bins, n_ctrl, cells, exclude` — no seed. `set.seed()` appears nowhere
in `R/malignancy.R`, and the k-means at line 703 uses `nstart`. Measured on one
fixed input, two calls differ by **0.3377** in module score. The formal named
`seed` (`R/malignancy.R:376`, documented at `:312`) is the output of
`FindMalignantSeedGenes()` — so a name-based seed audit reports this function as
seeded when it is the one unseeded path. The other eight seeded paths are correct
and do record their seed (`reduction$umap$params$seed`,
`clustering$louvain_clusters$params$seed`).

**3. Two tests assert on source text instead of behaviour — 10 assertions that
break on a renamed internal.** `tests/testthat/test-ai-evidence-rerun.R:70-72`
greps `deparse(body(fn))` for `sclet_ai_record_evidence_for_run` across 8
recorders; `tests/testthat/test-env-path-version-free.R:91-92` matches
`LD_LIBRARY_PATH=` inside `deparse(body(sclet_env_run))`. Both survive
`R CMD check` (they are location-independent, so this is *not* the old pin-test
trap) but they pin spelling, not behaviour. The behavioural replacement is already
written in the same file — `test-ai-evidence-rerun.R:41-44` expects
`evidence_id_conflict`, which is the observable consequence of routing through the
shared recorder. Loop that over the eight names instead of grepping bodies.

### What I would fix first, in order

1. Tests for the README quickstart path (`Read10X`, `subset_feature`,
   `PercentageFeatureSet`, `QCMetrics`) — cheapest, highest signal.
2. A `seed_rng =` argument threaded through `sclet_add_module_score()` and
   `ClassifyMalignantEpithelium()`, recorded in `params` like the other paths.
3. Rewrite the 10 source-text assertions as behaviour assertions.
4. `^Rplots\.pdf$` into `.Rbuildignore` (one line).

### What I could not fault

- DESCRIPTION: no undeclared package; no unused declared package (checked both
  directions — R CMD check only checks one).
- Clean-checkout build: every `R/` and `man/` file tracked; no untracked build
  input; installs and checks with a single deliberate NOTE.
- Test-file duplication: none, measured three ways.
- README: 67/67 named functions exported, all documented arguments real formals,
  quickstart runs.

### Not verified (and why)

- `devtools::test()` (15 min) and `make check` (25 min) were not run — the brief
  forbids the full check, and a tarball-level `R CMD check` answered the packaging
  question in 5 min.
- Python-backed stochastic paths (`RunCellRank`, `RunCellFate`, `RunVelocity`,
  `RunRegVelo`, `RunDecontX`, `RunSuperCell`) — need basilisk environments,
  outside the cheap tier; seed propagation into `scvelo`/`cellrank` is unvalidated.
- README step 5 (enrichment) — needs an annotation source; avoided a network call.
- The `.git`-in-tarball mechanism (§4.2-b) — seen once, not reproducible.
