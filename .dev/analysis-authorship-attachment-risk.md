# Analysis: the sole-node attachment rule is unsafe as an extension pattern

Date: 2026-10-13
Author: Lead (measured directly, not inherited from a subagent report)
Status: finding, blocks the direction of `.dev/spec-ai-evidence-authorship.md` section 6

## 1. What I tested, and what came out

`spec-ai-evidence-authorship.md` section 7 named one risk in the abstract:

> attaching to a finding whose *statement* is unrelated to the evidence (for
> example a finding about cell-cycle phase on an object whose only node is
> about markers). The one-node rule cannot detect that.

It is not abstract. Two offline measurements, both reproducible:

**(a) The rule attaches an unrelated finding, and raises it.**

A single `ev:findallmarkers` node (`kind = "test"`, `claim_level =
"associated"`) and two findings, both citing nothing:

| finding | subject | requested | recorded | rule |
|---|---|---|---|---|
| "marker counts describe recorded test output only" | the marker test | hypothesis | **associated** | `sole_evidence_node` |
| "no mitochondrial or ribosomal content check and no doublet detection were recorded" | the ABSENCE of QC records | hypothesis | **associated** | `sole_evidence_node` |

Audit reported `status = "downgraded"`, `n_problems = 0`. The second finding is
about evidence that does **not** exist; the engine supported it with a marker
test. That is the system manufacturing a basis for a claim, in the one
direction the apparatus exists to block.

**(b) Attaching more nodes LOWERS the ceiling, so "attach them all" is not an
extension.** The ceiling is the minimum rank over the cited nodes within each
dependency group (`R/ai-claims.R`):

| cited set | independent groups | ceiling |
|---|---|---|
| `ev:strong` alone (`consistent_with`) | 1 | `consistent_with` |
| `ev:strong` + `ev:weak` (`observed`, same group) | 1 | **`observed`** |
| `ev:strong` + independent `consistent_with` (distinct scope) | 2 | `consistent_with` |
| the above three together | 2 | **`associated`** |

Adding a node inside an already-covered group can only attenuate. Any
multi-node rule that attaches a *set* would therefore degrade claims relative to
attaching nothing, and relative to attaching the one relevant node.

**(c) How often (a) would bite, bounded.** Across the four
`interpret_evidence_refs.jsonl` records, 51 of 60 findings cite nothing and
would all be attached by the sole-node rule on a one-node object. **46 of 51
(90%) never name the node anywhere in their own text.**

The 90% is an **upper bound on the false-attachment rate, not a measurement of
it.** A finding can rest on evidence without typing its id (the whole reason
authorship exists is that models under-cite). What is measured is that the rule
has no relevance check at all, and that the example in (a) — a finding about
missing QC records — is present in the real model's output, not invented.

## 2. Why this changes the extension, not just adds a caveat

Section 6 of the authorship spec anticipated the next step as widening the
mapping rule: *extend coverage with a defensible rule (for example matching a
finding's text against a node's `kind`, or using the `source` step's declared
outputs)*. Both candidate mechanisms fail against (a):

- **Text matching against `kind`** is keyword matching on prose. The finding in
  (a) row 2 mentions "doublet detection"; a `test`-kind node from a marker run
  would match nothing or match wrongly. Nothing about `kind` tells you which
  claim a sentence is making.
- **`source` step output matching** requires knowing which step a finding
  discusses — which is the same unknown, moved one level down.

So neither yields a *deterministic* mapping, and determinism is what makes
authorship defensible at all: the engine may only assert what it can compute.

## 3. What the engine actually knows -- and the fix this section first proposed, falsified

**Correction (measured after the first draft of this document).** I proposed
that the defensible extension was *scope rather than select*: attach only where
the interpretation is confined to one execution, because then the question and
the evidence share a boundary and the engine need not read the finding's mind.
That proposal does not survive a direct test. Evidence-node visibility by
`target` on a two-node object (`ev:findallmarkers` from a marker test,
`ev:traj` from a summary):

| `target` | nodes visible to the model |
|---|---|
| `NULL` / `report` / `status` / `qc` | both (whole-object review) |
| `markers` | `ev:findallmarkers` only |
| `trajectory` | none |
| an execution id, e.g. `plan_analysis_20261013_...` | **none** |
| `run_de_test` (the producing action's own name) | **none** |

`AIExplainAnalysis(object, id = execution_id)` passes `target = id`
(`R/ai-tools.R:151`), and `sclet_ai_target_matches()` filters evidence by
matching `paste(type, id)` against alias/keyword patterns. An execution id
matches no node, because nodes are tagged with the *action* that produced them
(`source`), not with the execution that ran it. So the scoped interpretation
path -- `RunAIAnalysis(interpret = TRUE)` -- sees an **empty evidence set**: it
is not a place where attachment could be made safely, it is a place where
nothing can be attached at all. Scoping by category (`markers`) works, but the
interpretive entry points do not use it. Category scoping is therefore not
available as a relevance mechanism on the paths that need it; making it so is a
ledger/provenance change (tag nodes with the execution id that produced them),
not an attachment-rule change.

What remains true from the earlier draft:

- Each node records `source` -- the action name that produced it
  (`findallmarkers`-style), validated against a completed analysis at
  `RecordAIEvidence()` time (`sclet_ai_evidence_source()`).
- `ExecuteAIPlan()` keeps per-step `results[[i]]$id/$action` and each step's
  `output_summary`, so nodes created during an execution *could* be
  attributed, if they were tagged.

## 4. Consequences to carry into the next spec

1. **Do not widen the sole-node rule to a most-relevant-node rule.** (b) says a
   wrong extra attachment is worse than none, and (a) says "wrong" is reachable
   today. Relevance cannot be computed from `kind` or `source` alone.
2. **The attachment is wrong on the one path that ships it.** `qc` and `status`
   are whole-object reviews: every node is in scope for every finding, so there
   is no confined subject. On the measured real data, 3 of 7 findings in a
   `qc_review` run and 10 of 16 in a `status` run cite nothing, and the
   attachable ones are dominated by findings about **absent** records -- the
   case where attaching an existing unrelated node is provably wrong:

   | rec | finding | what it claims | cites | names the node |
   |---|---|---|---|---|
   | 3 | R3 | "no recorded mitochondrial/ribosomal content check, no doublet detection, no cell-level filtering summary" | no | no |
   | 3 | R6 | "mainline_missing = [velocity_inputs, velocity, fate, perturbation]" | no | no |
   | 3 | R7 | "blocked_actions.execute_analysis.blocked = null ... available_actions is empty" | no | no |
   | 4 | f1 | "Object dimensions and available assays/reductions/graphs" | no | no |
   | 4 | f10 | "Missing prerequisites reported by the status payload. Also absent: trajectory..." | no | no |

   R3/R6/R7/f10 assert that evidence is *missing*. The current rule hands them
   `ev:findallmarkers` and raises `observed` -> `associated`. The audit reports
   this as a clean `downgraded` correction with `n_problems = 0`, so nothing
   downstream can tell it apart from a legitimate raise.
3. **Therefore the shipped rule needs a replacement, not an extension.** Three
   candidates, measured:
   - (i) **Never attach to a finding whose text claims something is missing.**
     A keyword guard. Brittle, and negative-only: its failure mode is the old
     floor rather than a wrong raise, which is why it is at least safe -- but it
     still attaches `ev:findallmarkers` to f1 ("Object dimensions and available
     assays/reductions"), which is merely less wrong.
   - (ii) **Stop attaching on whole-object review targets** (`qc`, `status`).
     Given the falsification in section 3 that means *nowhere* today: a revert
     of authorship until node provenance is tagged with the execution that made
     it. Correct-by-construction, and it forfeits the upgrades in section 9.
   - (iii) **Complete a citation the model already made.** Attach only when the
     finding's own prose names the node's id while `evidence_refs` is empty.
     Measured on the same four records: **5 of the 51 attachable findings
     (10%)** do exactly that -- e.g. a `status` finding reading
     "AI evidence ev:findallmarkers: status completed, method ..." with
     `evidence_refs=` empty. The model identified the evidence and only failed
     to fill the field; the engine moves the id from prose to structure. It is
     deterministic, it cannot fire on R3/R6/R7/f10 (those never name the node,
     they name what is *missing*), and it needs no provenance change.

   **Recommendation: (iii) now, and (ii)'s discipline kept as the rule for
   anything wider.** (iii) is a 10% yield against (a)'s 90% -- the honest
   framing is that authorship mostly does not work, and the working part is
   transcription, not attribution. (ii) should still be the fallback if (iii)
   cannot be landed cleanly. Provenance tagging (nodes carry the execution id
   that produced them, so an execution-scoped interpretation can see its own
   evidence) is the real enabler for any future widening and belongs in its own
   spec.
4. **The online batch cannot detect this at all.** Every record in it has one
   node, so `ambiguous_multiple_candidates` never fired; and `audit_status =
   downgraded` looks healthy whether the attachment was relevant or not. The
   numbers written into `spec-ai-evidence-authorship.md` section 9 (6/7, 12/13,
   9/10 upgrades) therefore **include** this defect firing -- they are not clean
   evidence that authorship works, and the spec must say so. A multi-node
   benchmark (`authorship_multinode.R`) must build the *irrelevant-finding* case
   in, or it will report success on exactly the behaviour this analysis rejects.
