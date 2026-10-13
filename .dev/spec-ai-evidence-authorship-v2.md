# Spec v2: attach only what the finding names, and clamp instead of overwrite

Status: ready for implementation
Baseline: `cd1427a` (all `test-ai-*.R` pass; `make check` 0/0/0)
Supersedes: the attachment rule and the derivation rule of
`.dev/spec-ai-evidence-authorship.md` (v1) and `.dev/spec-ai-derived-claim-levels.md` §4a step 4
Evidence: `.dev/evidence-node-census.md`, `.dev/analysis-authorship-attachment-risk.md`,
`.dev/analysis-claim-raise-vs-clamp.md`

## 1. What v1 got wrong, measured

Three statements in v1's own spec text are contradicted by measurement. They are
kept here because they are what a future reader needs to not re-assume.

1. **"A real analysis records several evidence nodes"** (v1 §6). False. The
   census and all 17 real-provider records that carry the count agree: the
   eight pipeline actions record **zero** nodes, and a standard
   normalize -> HVG -> scale -> PCA -> [UMAP] -> neighbors -> clusters ->
   `run_de_test` run ends with **exactly one** (`ev:findallmarkers`) in 7 of 7
   executed records. So v1's refusal-to-guess safety property almost never
   fires, and its attach branch -- the unsafe one -- is the common case.
2. **v1 §4 req. 2 required a relevance check the code never implemented.** The
   spec said the node's `source` must name "a step the finding's text does not
   contradict"; `sclet_ai_attach_finding_evidence()` filters on completed /
   non-decision / exactly-one and never reads the finding's text (verified by
   reading the function). The consequence is measured: on a one-node object, a
   finding asserting that evidence is **absent** ("no recorded
   mitochondrial/ribosomal content check, no doublet detection") is attached
   `ev:findallmarkers` and raised `hypothesis -> associated`, with audit
   `n_problems = 0`. Across the four real records, 46 of 51 cite-nothing
   findings (90%) never name the node anywhere in their prose; the upper bound
   is stated as such in the risk analysis because not-naming is not the same as
   not-resting-on.
3. **The fix this session first proposed (scope attachment to the producing
   execution) is falsified.** Evidence visibility is keyed by `target` matched
   against `paste(type, id)`; `target = "markers"` shows the marker node, but
   `target = <execution id>` or the producing action name shows **no** nodes at
   all, because nodes are tagged with the action `source`, never with the
   execution. So the execution-scoped interpretation path is an empty evidence
   set, not a safe place to attach.

Separately, the derivation's raise rule (v1's predecessor, §4a step 4) converts
a model's `observed` self-report into `associated` -- a rank-1 "I looked at this"
into a rank-3 relational claim -- and did so 15 times in the 2026-10-13 batch,
including `measured -> associated` x3. The audit reports every such move as
`entry_status = "downgraded"` and counts it in `n_downgraded`, so an upgrade is
indistinguishable from a downgrade in the logs. No existing test covers the
`observed`-request case, which is why it shipped.

## 2. The rule

**Attachment: transcribe, never attribute.** A finding's own prose already
identifies the evidence in the cases the engine can act on with certainty. So:

- model supplied `evidence_refs` -> unchanged (as v1);
- requested level `causal` or `estimated` -> never gain a ref (as v1);
- otherwise scan the finding's `statement` and `description` for exact
  evidence-node ids among the completed non-`user_decision` nodes:
  - **exactly one** id found -> attach that id (`prose_named_unique`);
  - none -> decline (`no_prose_reference`);
  - two or more -> decline (`ambiguous_prose_reference`).

This replaces v1's sole-node rule. It is deterministic (string equality against
ids the ledger owns, not keyword heuristics about meaning), it is the
**multi-node** extension v1 deferred -- selection comes from the model's own
words, so two nodes no longer force a decline -- and it deletes the C5 failure
class rather than narrowing it: R3/R6/R7 name no node, so they keep the floor.
Measured yield on the existing records: 5 of 51 cite-nothing findings qualify
(10%), against 51 under v1's rule. The lower yield is the point; the other 46
were the defect firing.

**Derivation: clamp, and raise only from the floor.** Replace "derived level is
the ceiling" with:

```
mapped   <- sclet_ai_finding_allowed_level(requested)      # suggestive -> consistent_with
req_rank <- rank(mapped), NA -> 0                         # hypothesis/absent = 0
final    <- if (req_rank == 0) ceiling_rank else min(req_rank, ceiling_rank)
```

so a request at the floor takes the ceiling (the raise is the only thing
authorship exists for, and under the new attachment rule it fires precisely on
findings that named their evidence), while a request already above the floor is
respected up to the ceiling. The alternatives were enumerated in the memo; two
are worth recording as rejected *with their measured cost*: `min()` with no
floor exemption makes attachment a pure no-op on the only objects it fires on
(it reproduces the 0-of-23 / 0-of-13 result), and ceiling-only (v1's rule)
manufactures the over-claim.

**Reporting: say which way it moved.** `entry$status` becomes `ok` | `raised` |
`corrected` (downward) | `overclaimed`, with `n_raised` and `n_lowered` counts,
so `audit_status` in the logs no longer calls an upgrade a downgrade. `causal`
keeps `overclaimed` + refusal; nothing about the ban or `RecordAIResult()`'s
refusal path changes.

## 3. Requirements

1. `sclet_ai_attach_finding_evidence()` implements §2: prose-name scan over
   completed non-`user_decision` nodes, unique-hit semantics, rule names
   `prose_named_unique` / `no_prose_reference` / `ambiguous_prose_reference` /
   `model_cited` / `requested_level_excluded` / `none`. The sole-node rule is
   removed, not kept as a fallback.
2. Node ids are matched as whole tokens against the prose, so `ev:markers` does
   not match inside `ev:markers_extra`. Keep the match case-sensitive (ids are
   machine-generated).
3. `sclet_ai_derive_finding_claim_level()` implements the clamp with the floor
   exemption. `estimated` must still return at step 1, before refs and rank are
   consulted (moving it after would rank 0 and raise a provenance tag into a
   strength claim).
4. `sclet_ai_audit_result_claims()` reports direction as in §2. The refusal
   semantics for `causal` are unchanged, and an attached-and-raised finding is
   `raised`, never `overclaimed`.
5. `evidence_refs_attached` / `attachment_rule` / `evidence_attachment_rule`
   keep their v1 meaning and stay persisted on the recorded finding; the new
   rule names simply replace the old string values.
6. `sclet_ai_findings_needing_support` (dead constant, referenced nowhere
   outside its own definition) is deleted.
7. No change to `sclet_ai_claim_ceiling()` (its min-over-nodes and
   independent-group arithmetic is what makes both rules safe), to
   `RecordAIEvidence()`, to the ledger, to any prompt, or to the `causal` ban.

## 4. Acceptance criteria

1. The C5 case: one node `ev:findallmarkers` (`associated`) + a cite-nothing
   `hypothesis` finding whose text asserts absence of records -> **not
   attached**, recorded `hypothesis`, `attachment_rule = no_prose_reference`.
2. The transcription case: the same object + a cite-nothing finding whose prose
   contains `ev:findallmarkers` -> attached, recorded `associated`
   (raise-from-floor), `attachment_rule = prose_named_unique`.
3. Multi-node selection: two nodes, prose names exactly one -> that one is
   attached, and the ceiling is computed from **that** node only (not from the
   pair, which would attenuate it: a second node in the same dependency group
   lowers the ceiling -- measured in the risk analysis).
4. Multi-node ambiguity: two nodes, prose names both -> declined, floor holds.
5. Clamping: `observed` request + `associated` ceiling -> recorded `observed`,
   status `ok`, **zero** movement counts. `measured` + `associated` ->
   `measured`. `consistent_with` + `associated` -> `associated`, `corrected`.
6. Floor raise survives: `hypothesis` (or absent level) + a model-cited ref at
   `associated` -> `associated`, status `raised`.
7. `estimated` is untouched by both rules and gains no ref; `causal` is still
   refused, gains no ref, and still blocks recording.
8. The v1 sole-node rule no longer exists: a cite-nothing finding that does not
   name the node is left at the floor even when it is the only node on the
   object (this is criterion 1, asserted as a regression guard against
   re-adding the fallback).
9. `make check` 0/0/0, full suite passes, and every test updated to the new
   semantics names the reason in a comment (several v1 tests assert exactly the
   behaviour this spec removes).

## 5. What this does not do

- **No relevance reasoning.** It cannot tell whether a naming is sincere; it
  only refuses to guess. A finding that names `ev:findallmarkers` while
  describing something else would still be attached -- but that requires the
  model to have identified the node, which is the same evidence a human reader
  would use.
- **No provenance tagging.** Nodes still do not carry the execution id that
  produced them, so no path can confine evidence to one execution. That is a
  separate, larger change (§1 point 3) and must get its own spec.
- **No new measurement of model behaviour.** The 10% yield in §2 is computed
  from four existing records on one-node objects; multi-node prose-naming is
  unmeasured because no real record has ever had two nodes. The benchmark
  action must produce a >= 2-node object or state plainly that it could not.
- **No prompt change.** As v1 §3: teaching the model to cite in prose is a
  different intervention with its own risk, and the structured `evidence_refs`
  channel from `3cd9af6` remains the primary route.

## 6. Risk

The rule's safety property is now "the model named it", which is verifiable but
also **quiet**: if the model stops naming ids in prose, attachment stops
entirely and the ledger goes back to `hypothesis` everywhere, indistinguishable
from a regression in the derivation. Mitigation is reporting: `n_raised` and the
`attachment_rule` counts make "0 attempts" visible, and the benchmark records
both. The opposite risk -- over-attaching -- is what v1 had, and it is measured.
