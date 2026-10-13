# Spec: engine-attached evidence authorship

Status: implemented and measured (2026-10-13)
Baseline: `1bddedc` (all `test-ai-*.R` pass; `make check` 0/0/0)

## 1. The open question this closes

`.dev/spec-ai-derived-claim-levels.md` built the derivation and measured it as
necessary but not sufficient (§4c). The derivation is correct and tested: given
`evidence_refs`, `sclet_ai_derive_finding_claim_level()` consults the ceiling
and raises the finding to the level the cited evidence can carry. It found zero
upgrades because **the model supplied no refs**, and the derivation's "no refs"
path correctly lands every finding at `hypothesis`.

§4c left two routes and rejected a third:

> do not instruct the model to cite evidence and claim more strongly, because
> its current non-fabrication is a measured virtue and that instruction trades
> it for a number. The plausible routes are (a) have the engine attach the
> evidence a finding relates to, since it knows which evidence each step
> created, leaving the model only to describe; or (b) accept that `claim_level`
> stays at the floor for interpretive reviews and stop treating the field as
> informative.

This spec implements **(a)**.

## 2. Measurement, refreshed

The §4c measurement predates the `3cd9af6` schema fix that gave `evidence_refs`
an output channel. Re-read directly from
`sclet-benchmark-ai/results/derived_claim_levels.jsonl` (four records,
2026-10-13):

| condition | task | requested | recorded | up | down | ceiling |
|---|---|---|---|---|---|---|
| `with_evidence` | `qc_review` | hypothesis x23 | hypothesis x23 | 0 | 0 | `allowed` |
| `with_evidence` | `status` | hypothesis x13 | hypothesis x13 | 0 | 0 | `allowed` |
| `without_evidence` | `qc_review` | hypothesis x17 | hypothesis x17 | 0 | 0 | `unsupported` |
| `without_evidence` | `status` | hypothesis x17 | hypothesis x17 | 0 | 0 | `unsupported` |

The two `with_evidence` rows are the ones that matter: `derivation_available`
is `TRUE`, the object carried `ev:findallmarkers`, the ceiling was `allowed`,
and still nothing moved — because the finding carried no `evidence_refs` and the
"no refs" path is the floor. So the missing ingredient is precisely the
reference, not the derivation. §4c's diagnosis holds and this is still the right
thing to build.

A later measurement (`interpret_evidence_refs.jsonl`) shows the schema fix did
make the model cite **sometimes**: 3 of 7 findings in one `qc_review` run and 6
of 16 in one `status` run. So the model can cite; it cites unreliably. That is
the gap authorship closes.

## 3. The design problem, stated honestly

Attaching evidence on the model's behalf is a **new way to be wrong**, and it is
more dangerous than the problem it solves. An over-claiming AI finding is the
failure this whole apparatus exists to prevent, and attaching evidence raises a
claim (`hypothesis` -> the evidence's own level). If the engine attaches
evidence that does not actually support the finding, the audit will report a
clean result for an unsupported claim — strictly worse than today's accurate
"everything is a hypothesis".

Two properties of the existing machinery bound the damage, and the design leans
on both:

1. **The ceiling is the minimum over the cited nodes' own claim levels**, and it
   reasons about independent lines of support (`n_independent_groups`). Taking
   it can never assert more than the evidence already declares for itself. So
   authorship cannot invent strength out of nothing — it can only surface
   strength that a recorded node already carries.
2. **Evidence nodes record their own provenance**: `source` (which completed
   analysis created them), `kind` (`deterministic_summary`, `test`, `plot`,
   `state`, `user_decision`), `parents` and `scope`. The engine therefore knows
   which step created each node, which is the basis of §4c's "it knows which
   evidence each step created".

Therefore the governing rule is: **attach only when the mapping from finding to
evidence is forced, and never guess.** A finding the engine cannot confidently
map stays at the floor. Under-attaching is the safe direction and is expected.

## 4. Requirements

1. Add `sclet_ai_attach_finding_evidence()` (internal). Given the object and one
   finding, it returns the finding with `evidence_refs` filled **only** when
   exactly one candidate is defensible, or the finding unchanged.
2. The rule for "forced" must be narrow and stated:
   - If the finding already carries `evidence_refs`, they are left exactly as
     given. The model's own citation is never overwritten.
   - Otherwise, if the object carries **exactly one** completed evidence node
     that is *not* a user decision, and that node's `source` names a step the
     finding's text does not contradict, the node is attached.
   - Otherwise nothing is attached.
3. **A one-node object is the only case this spec attaches on.** With two or
   more candidate nodes the mapping is a guess, and the correct answer is to
   leave the finding at the floor. This deliberately makes the feature a no-op
   on objects with several evidence nodes; that is a consequence of choosing
   under-attachment, not an oversight, and it must be documented as such.
4. `kind = "user_decision"` nodes are never auto-attached. The ceiling already
   separates them, and a finding is not "supported" by a decision the user made.
5. The attachment must be **recorded, not silent**: the audit entry for an
   attached finding must show it was attached by the engine rather than cited
   by the model. Two new fields on the derivation's return value:
   `evidence_refs_attached` (logical) and `attachment_rule` (character, the rule
   that fired, or `"none"`).
6. Authorship changes where a level comes from, never how strict the bar is:
   no change to `sclet_ai_claim_ceiling()`, `sclet_ai_derive_finding_claim_level()`'s
   rule (still not `min(requested, ceiling)`), `sclet_ai_audit_result_claims()`'s
   verdicts, or the `causal` ban.
7. `RecordAIResult()` is where this is applied, immediately before the existing
   derivation call, so a finding that gains a ref gains its derived level in the
   same pass. No new user-facing argument: authorship is not opt-in, because a
   correct attachment is strictly more informative than none and the floor is
   unchanged when it does not fire.
8. No prompt change. §3 of the earlier spec still stands.

## 5. Acceptance criteria

1. A finding with no refs, on an object with exactly one non-user-decision
   evidence node, gains that node's id; the returned derivation reports
   `evidence_refs_attached = TRUE` and a non-`"none"` `attachment_rule`, and the
   recorded level is the evidence's own level rather than `hypothesis`.
2. The same finding on an object with **two** evidence nodes is left at
   `hypothesis` with `evidence_refs_attached = FALSE`. This is the
   anti-guessing criterion and the most important one in this spec.
3. A finding that already carries `evidence_refs` keeps them verbatim, even when
   the object has exactly one node and it differs. The model's citation is
   authoritative for what it cited.
4. A `user_decision` node is never attached, even when it is the only node.
5. A finding requesting `causal` is still refused and never gains a ref — the
   ban is unaffected by authorship.
6. `estimated` findings keep their provenance level and gain no ref (the
   provenance path returns before evidence is consulted).
7. A finding attached by the engine is reported as attached in the audit entry,
   so a reader can tell authorship from a model citation.
8. With authorship disabled or inapplicable, every existing `test-ai-claims.R`
   assertion still holds — the change must not alter any verdict on a result
   that already cited its evidence.
9. `make check` 0/0/0.

## 6. What this deliberately does not do

- **No guessing with multiple evidence nodes.** The feature will therefore do
  nothing on the pbmc3k `with_evidence` measurement above, since that object
  carries one node — but a real analysis records several
  (`ev:qc`, `ev:clusters`, `ev:findallmarkers`). This spec accepts that, because
  the alternative is attaching the wrong one. Extending coverage with a
  defensible rule (for example matching a finding's text against a node's
  `kind`, or using the `source` step's declared outputs) is the natural next
  spec, on evidence from a run with several nodes.
- **No re-run of the benchmark.** The measurement is the refreshed reading of
  the existing logs plus the unit-level contract.
- **No prompt change** and no instruction to cite more.
- **No change to route (b)** — if later measurement shows authorship cannot fire
  in practice, §4c's route (b) remains the honest fallback: stop treating
  `claim_level` as informative for interpretive reviews.

## 7. Risk

The main risk is a subtle one: attaching to a finding whose *statement* is
unrelated to the evidence (for example a finding about cell-cycle phase on an
object whose only node is about markers). The one-node rule cannot detect that,
so criterion 2 is a guard on count, not on meaning. This is why requirement 5
exists — the attachment is visible in the audit rather than laundered into a
clean-looking result — and why the feature must be described to users as
"the engine attached the only available evidence", not "this finding is
supported".

## 8. Outcome

Implemented on `devel`. Two things were measured while implementing that
the design above did not predict, and both changed the code:

1. **Authorship must run after the ban and provenance checks, not before.**
   The first wiring attached before consulting the derivation, which handed
   `causal` findings exactly the reference `validate_sclet_ai_result()`
   demands for a causal claim -- the ban defeated by helpfulness, and `estimated`
   findings a reference their level ignores by construction. Both levels now
   return `requested_level_excluded` before the node lookup (test-ai-evidence-
   authorship.R scenario 4).
2. **The recorded finding carries two similar names and `$` reads the wrong
   one.** `evidence_refs` (character) and `evidence_refs_attached` (logical)
   coexist on the recorded finding, and R's `$` partial-matches: a finding that
   cited nothing read back `FALSE` as its evidence id, `as.character()` made it
   the literal `"FALSE"`, and every such result failed to record ("evidence ref
   not valid: FALSE"). This was caught only by the full suite (4 tests in
   `test-ai-claims.R` / `test-ai-derived-claim-levels.R`), not by the new file.
   Fixed by reading with `[[` under a name guard (`sclet_ai_finding_refs()` and
   the derivation's ref read); the regression guard asserts the honest
   `no_refs` verdict, not the phantom `unsupported` (scenario 12). Four older
   tests that pinned "sole node, no refs -> hypothesis floor" were split into
   the attachable case (now raised to `associated`) and a two-node case (floor
   holds), so both halves of the new contract stay pinned.

Acceptance criteria 1-9 all verified in
`tests/testthat/test-ai-evidence-authorship.R`: 12 scenarios, 52 assertions,
teeth confirmed twice (removing the attachment: 8 failures; reverting the
partial-match guard: 1 failure). Full suite 2805 pass / 0 fail / 78 warn
(unchanged) / 3 skip. `make check` 0 errors / 0 warnings / 0 notes.

The behaviour, measured offline on an SCE with one evidence node at
`associated` and a `hypothesis` finding citing nothing: audit entry
`claim_level = "associated"`, `evidence_refs_attached = TRUE`,
`attachment_rule = "sole_evidence_node"`, and the recorded ledger finding
carries both the raised level and the attached id. With two nodes:
`attachment_rule = "ambiguous_multiple_candidates"`, level stays
`hypothesis`, status `ok` -- the floor holds exactly as before.

No online re-measurement. Whether the real model's findings on pbmc3k are now
routed through the attached sole node (`interpret_evidence_refs` records show
`audit_status: downgraded` already, from the `3cd9af6` claim-level shifts) is a
separate billed question. Criterion 2 -- the anti-guessing behaviour -- is the
load-bearing one for that next measurement: on a real pipeline object with
several nodes, authorship must do nothing and report it.
