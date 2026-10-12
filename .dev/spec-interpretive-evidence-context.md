# Spec: let the interpretive context carry evidence ids

Status: accepted for implementation
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
   contains `ev:findallmarkers` as a string.
2. Same for `target = "status"`.
3. Same size-order of context (no large expansion beyond the evidence entries).
4. `sclet_ai_sanitize_context()` produces no new warnings and no new redactions
   on the `qc`/`status` contexts, compared with before.
5. `AIReviewQC()` and `AIStatus()` on the same object produce findings whose
   `evidence_refs` can be non-empty -- i.e. the vocabulary is at least available
   to the model.
6. `make check` 0/0/0, and the existing `test-ai-privacy.R` and
   `test-ai-claims.R` still pass.

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
