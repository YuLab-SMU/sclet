# Spec: let the interpretive context carry evidence ids

Status: implemented and measured (2026-10-13)
Baseline: `9c24bd0`

## 1. Problem, measured

The claim-ladder chain is: a finding cites `evidence_refs` -> the engine derives
the ceiling -> the recorded level is a function of the evidence. Measured
2026-10-12: **the model cites `evidence_refs` in 0 of 12 and 0 of 19 findings**
(`sclet-benchmark-ai/results/derived_claim_levels.jsonl`), so every finding
lands at `hypothesis` and `claim_level` stays a constant.

Before proposing a fix at the model level, the audit of 2026-10-12 required
checking whether the model can even see what it is being asked to cite.
**It cannot.** Verified against a pbmc3k object carrying `ev:findallmarkers`:

| context target | size | mentions `ev:findallmarkers` |
|---|---|---|
| `status` (`AIStatus`) | 2,671 chars | **no** |
| `qc` (`AIReviewQC`) | 3,626 chars | **no** |
| `report` | 10,385 chars | yes |

The two tasks that produce interpretive findings are given no evidence id at
all, so "cite evidence" is not merely discouraged -- it is **impossible**. This
is the same failure class as the `success_criteria` chain: machinery exists,
and the model was never given the vocabulary.

## 2. Root cause, in code

`sclet_ai_target_matches(type, id, target)` in `R/ai-context.R` decides which
ledger/state records a given context target may see.

- `target = "report"` returns `TRUE` unconditionally -- sees everything.
- `target = "qc"` has aliases `qc, preprocess, clean, doublet, phenotype`.
  Neither `ai_evidence`, nor `clustering`, nor `detest` matches, so the evidence
  record is filtered out.
- `target = "status"` has an **empty alias vector**, so nothing matches and the
  view is empty for that group.

`AIReviewQC()` uses `target = "qc"`; `AIStatus()` uses `target = "status"`;
`AIExplainAnalysis()` passes through whatever `type`/`id` it is given (so it can
see evidence when the caller names it, and otherwise does not).

## 3. Requirements

1. `target = "qc"` must include `ai_evidence` records.
2. `target = "status"` must include `ai_evidence` records.
3. **Privacy unchanged.** The payload classifier and `include_artifacts` gating
   must behave exactly as before for every field that is not an evidence id.
   Evidence *values* (`summary$values`) must not be broadened onto the
   interpretive path -- only ids, kinds, claim levels and status, which is what
   a finding needs in order to cite.
4. Existing behaviour for `planning` (which currently has no evidence exposure
   either) is out of scope for this spec: the plan path has its own measured
   gap and its own change.

## 4. Acceptance criteria

1. `sclet_ai_context(sce_with_evidence, target = "qc", detail = "full")`
   contains `ev:findallmarkers` as a string. — **met**: the benchmark's
   `evidence_id_in_context` reads FALSE before, TRUE after (commit `174f3b1`),
   for both `qc_review` and `status`.
2. Same for `target = "status"`. — **met**: same run, `status`
   `evidence_id_in_context` FALSE -> TRUE.
3. Same size-order of context (no large expansion beyond the evidence entries).
   — **met**: the context grew by the evidence entries only; the bounded
   summary is what `174f3b1` restored and then re-scoped, and no target changed
   size class.
4. `sclet_ai_sanitize_context()` produces no new warnings and no new redactions
   on the `qc`/`status` contexts, compared with before. — **met by tests**: the
   privacy suite passes unchanged (`test-ai-privacy.R`,
   `test-ai-context-evidence.R`); `make check` 0/0/0.
5. `AIReviewQC()` and `AIStatus()` on the same object produce findings whose
   `evidence_refs` can be non-empty -- i.e. the vocabulary is at least available
   to the model. — **met, and this is the criterion that needed the second
   commit**: the first commit made the vocabulary *available* and the model
   still cited 0 of 14 and 0 of 23 findings, so `3cd9af6` had to name
   `evidence_refs` in the finding schema and instruct the model to fill it.
   After both: 3 of 7 (`qc_review`) and 6 of 16 (`status`) findings cite
   `ev:findallmarkers`, zero invented.
6. `make check` 0/0/0, and the existing `test-ai-privacy.R` and
   `test-ai-claims.R` still pass. — **met**: `make check` 0 errors / 0 warnings
   / 0 notes, full suite 2617 pass / 0 fail, including `test-ai-privacy.R`,
   `test-ai-claims.R`, `test-ai-evidence.R` and the new
   `test-ai-findings-schema.R` (5 scenarios, 36 assertions).

## 5. Risk

Exposing evidence ids to more targets widens the context. The evidence record
already passes through `sclet_ai_normalize_record()`, which the privacy layer
screens; the risk is therefore limited to *more targets seeing the same record*,
not new content. The acceptance criteria include a before/after sanitizer
comparison so this is measured rather than assumed.

## 6. Measurement plan

Benchmark: `R/actions/interpret_evidence_refs.R` in `sclet-benchmark-ai`,
measuring before/after whether `AIReviewQC`/`AIStatus` findings carry
`evidence_refs`, on the same pbmc3k object. The current baseline is 0/12 and
0/19.

## 7. Outcome (measured 2026-10-13)

Benchmark `sclet-benchmark-ai`, action `interpret_evidence_refs`, fixture
pbmc3k = 2700 cells x 32738 genes, real pipeline normalize -> HVG -> PCA -> kNN
-> louvain 7 clusters -> `FindAllMarkers` producing evidence node
`ev:findallmarkers`, model `deepseek-chat`, records in
`results/interpret_evidence_refs.jsonl`.

| run | evidence id in context | findings citing `ev:findallmarkers` | reason |
|---|---|---|---|
| before (motivating run, earlier 2026-10-13) | FALSE | 0 of 14 (`qc_review`), 0 of 23 (`status`) | -- |
| intermediate (`174f3b1` only) | TRUE | 0 of 14, 0 of 23 | `model_did_not_cite` |
| after (`174f3b1` + `3cd9af6`) | TRUE | 3 of 7, 6 of 16 | `model_cited_correctly` |

What each step actually proves:

- The **before** run is the spec's premise reproduced on a live pipeline: with
  no evidence id in context, citation is impossible, not merely discouraged.
- The **intermediate** run is the important negative result. Exposing the ids
  was **necessary but not sufficient**. The model described the record
  correctly in prose -- "the only evidence record (ev:findallmarkers) is
  claim_level = 'associated', kind = 'deterministic_summary'" -- and still
  cited 0 of 14 and 0 of 23 findings. The cause is mechanical:
  `evidence_refs` was declared `z_array(z_any())` in
  `sclet_ai_native_result_schema()`, so the field had no output channel and the
  model had nowhere to put the id it clearly knew.
- The **after** run closes it. `3cd9af6` names `id`, `statement`,
  `description`, `severity`, `claim_level`, `evidence_refs` in the finding
  schema (enums shared with the local validator via
  `sclet_ai_finding_claim_levels()` / `sclet_ai_finding_severities()`), and
  `sclet_ai_evidence_citation_instruction()` tells the model to fill it with
  ids present in the supplied ledger, verbatim. `citation_improved: TRUE` and
  `reason: model_cited_correctly` for both targets.

Zero invalid citations: every cited id resolves to the real evidence node. The
citations are precise rather than blanket. Findings that are structural truths
-- dimensions, preprocessing methods, missing prerequisites -- carry no
`evidence_refs` at all. Every finding that rests on the FindAllMarkers summary
cites the node, and some state its limit in the finding text itself ("Claim
level is associated, not causal"; "No causal, trajectory, velocity, or
perturbation conclusions are supported by the recorded evidence"). Claim levels
stop being a constant: the before run put all 14 and all 23 findings at
`hypothesis`, the after run splits them into observed / associated / hypothesis
(qc: 5 observed + 1 associated + 1 hypothesis; status: 13 observed + 3
associated), with `audit_status: downgraded` where the model over-reached. The
claim-ladder chain the spec set out to repair now has live inputs.

One collateral finding, worth recording because it was not in the plan: the
in-progress context change had bounded `ai_evidence` summaries to
id/kind/claim_level/status for any non-`report` target, which stripped
`$source`, `$values`, `$dependency_group` and `$scope` from *internal* ledger
readers. `RecordAIEvidence`, `ValidateAIEvidenceRefs` and
`sclet_ai_claim_ceiling` all call `GetAnalysisLedger()` with `target = NULL`,
so ~222 tests failed. `sclet_ai_target_exposes_full_evidence()` in
`R/ai-context.R` returns TRUE when no target is requested, restoring full
evidence for internal consumers while the bounded view stays in force on the
interpretive path.
