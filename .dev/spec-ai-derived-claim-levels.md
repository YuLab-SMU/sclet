# Spec: derive finding claim levels from evidence instead of asking the model

Status: implemented (`a0a4c82`), measured, and found necessary but NOT sufficient -- see section 4c
Baseline: `4051f41` (sclet clean; all `ai-claims` tests pass)

## 1. Problem, measured

`sclet` has a claim-calibration apparatus: each AI finding carries a
`claim_level`, and `sclet_ai_audit_result_claims()` grades it against a
deterministic evidence ceiling computed by `sclet_ai_claim_ceiling()`. The
apparatus is implemented, unit-tested (18 tests in `test-ai-claims.R`), and
reached from `RecordAIResult(audit_claims = TRUE)`.

**The model never engages it.** Measured 2026-10-12 across seven real runs and
**114 findings** (`sclet-benchmark-ai`, `results/interpret_claims.jsonl`,
`results/interpret_null_signal.jsonl`):

| fixture | evidence | findings | claim levels |
|---|---|---|---|
| null (no signal) | none | 12 + 14 | hypothesis x26 |
| pbmc3k | none | 19 + 20 | hypothesis x39 |
| pbmc3k via `ExecuteAIPlan` | `ev:findallmarkers` | 10 + 27 + 12 | hypothesis x49 |

Not one finding above the floor. `hypothesis` is defined in `R/ai-claims.R` as
"explicitly speculative and self-limiting ... it asserts nothing", and
`sclet_ai_findings_needing_support` is `consistent_with`, `suggestive`,
`causal` -- so a hypothesis needs no evidence and **can never be graded as an
over-claim**.

Two consequences:

1. The clean audit results are *trivially* clean. They show the model does not
   fabricate, which is genuinely reassuring, but they do **not** show
   calibration, because nothing was ever at stake.
2. **`claim_level` carries no information today -- it is a constant.** Every
   consumer of it (`validate_sclet_ai_result()`'s rule that `causal` requires
   `evidence_refs`; the ceiling audit; any downstream report logic) is testing a
   field that never varies.

This is the same failure shape as the `success_criteria` chain: well-tested
machinery the model never reaches for. The difference is that there the wiring
was broken at four points and fixing the wiring worked. Here the wiring is
intact, and **asking the model harder is the wrong lever**: the model is
currently not fabricating, and the most likely effect of instructing it to use
stronger levels is to induce over-claiming -- converting a real strength into a
risk.

## 2. Proposal

Stop asking the model to self-report a claim level. **Derive it from evidence,
deterministically, and keep the model's level only as a request that the engine
clamps.**

The engine already computes the answer and throws it away. Verified on a pbmc3k
object carrying `ev:findallmarkers`:

```
sclet_ai_claim_ceiling(sce)
  status: allowed   ceiling: associated   n_independent_groups: 1
  reasons: all supporting evidence resolves to 1 independent line(s) of
           support; a claim cannot exceed 'associated'

sclet_ai_claim_ceiling(sce, proposed_claim_level = "consistent_with")
  status: downgraded   ceiling: associated

sclet_ai_claim_ceiling(sce, proposed_claim_level = "observed")
  status: allowed      ceiling: associated
```

and, end to end through the audit today:

```
finding claim_level = "consistent_with", evidence_refs = ["ev:findallmarkers"]
  -> status: "overclaimed", n_problems: 1, allowed_claim_level: "associated"
```

Read that last line carefully, because it is the whole argument. The audit has
**already computed the correct claim level** (`allowed_claim_level =
"associated"`) and still reports the result as a failure. The engine knows what
the finding should say and declines to say it. This is not a request to build a
derivation -- the derivation exists and is tested -- it is a request to *use* it
instead of only complaining about the absence of it.

So `claim_level` becomes an **output of the engine**, not an input from the
model, and it varies with the evidence instead of being constant.

### 2.1 The vocabulary mismatch, which is the real design problem

The two ladders are different, and the spec must not paper over it:

| | levels |
|---|---|
| evidence (`sclet_ai_claim_levels`) | `observed`, `measured`, `associated`, `consistent_with` |
| finding (`validate_sclet_ai_result`) | `observed`, `measured`, `estimated`, `associated`, `consistent_with`, `suggestive`, `hypothesis`, `causal` |

The ceiling can only ever return one of the four evidence levels. The finding
ladder's extra levels mean different things:

- `estimated` -- a *provenance* statement (a value was modelled, not measured),
  not a strength statement. It is orthogonal to the ceiling and must survive it.
- `hypothesis` -- the floor, asserting nothing. A finding that the engine cannot
  support at all should land here, not at `observed`.
- `suggestive` -- one step above `consistent_with` on the finding ladder, and
  mapped *down* to `consistent_with` by `sclet_ai_finding_allowed_level()`
  before the ceiling is consulted. Already handled.
- `causal` -- banned outright for AI findings. Unaffected by this proposal.

So the derivation rule needs to be explicit about where each level comes from.

### 2.2 Derivation rule

For each finding, with its `evidence_refs` (possibly empty):

1. If the model's level is `estimated`, keep it -- provenance, not strength.
2. If the model's level is `causal`, it stays a problem (`causal` is never
   assertable). Do not silently downgrade it; the audit must still say so.
3. Compute `ceiling <- sclet_ai_claim_ceiling(object, evidence_ids = refs)`.
   - **No refs, or `ceiling$status == "unsupported"`** -> derived level is
     `hypothesis` (the floor that asserts nothing). This is the correct landing
     place for an unsupported finding, and it is strictly safer than today's
     behaviour of leaving the model's word in place. Note this is *not* what the
     audit currently suggests: with no refs it reports
     `allowed_claim_level = "observed"`, which would over-state an unevidenced
     finding. `hypothesis` is the right floor.
   - **Otherwise** -> derived level is
     `min(model's level mapped through sclet_ai_finding_allowed_level(),
     ceiling$ceiling_claim_level)` by `sclet_ai_claim_rank()`.
4. Record both: `claim_level` (the derived value) and `requested_claim_level`
   (what the model asked for), so a downgrade is visible rather than silent --
   the same honesty rule the ceiling function already follows with its
   `status = "downgraded"`.

Reported in the audit as a `downgraded` entry when they differ, so
`RecordAIResult(audit_claims = TRUE)` can surface it without refusing the result
(a downgrade is a correction, not a violation).

### 2.3 Where this gets applied

The derivation belongs in the same place as the existing audit, so there is one
definition:

- `sclet_ai_audit_result_claims(object, result)` gains the derived level in its
  per-finding `entry` (`claim_level`, `requested_claim_level`,
  `derived_from_ceiling` logical).
- `RecordAIResult(object, result, audit_claims = TRUE)` writes the **derived**
  level into the finding it persists, never the requested one, so the recorded
  level is the derived one.

  Verified accessor note: the findings reach the ledger via
  `sclet_get_state(object)$analyses[[id]]$artifacts$findings`, **not** via
  `GetAnalysisLedger(object)$analyses[[id]]` -- the latter strips `artifacts`
  unless `include_artifacts = TRUE`, so a naive read-back shows no findings at
  all and the measurement would silently observe nothing. Anything verifying a
  recorded claim level must use the state accessor (or pass
  `include_artifacts = TRUE`).
- `validate_sclet_ai_result()` keeps its shape-only checks (it has no object, so
  it cannot derive anything) and continues to reject an invalid level *name*.

## 3. What this is not

- **Not a prompt change.** Adding "use stronger claim levels when you have
  evidence" to the system prompts is the obvious alternative and is rejected
  above: the model's current non-fabrication is a property worth keeping, and
  that instruction trades it for a number that the engine can compute itself.
- **Not a relaxation of the audit.** `causal` stays banned, unsupported levels
  stay reportable. The change is where the level comes from, not how strict the
  bar is.
- **Not implemented here.** This spec is deliberately left at proposal status.

## 4. Acceptance criteria

Each of these was checked against the current code before being written, so the
spec states observed behaviour rather than assumed behaviour. "Today" is what
`sclet_ai_audit_result_claims()` does now; "after" is the required behaviour.

| # | scenario | today (verified) | after |
|---|---|---|---|
| 1 | no `claim_level` | `hypothesis`, status `ok` | unchanged |
| 2 | `estimated` | passes untouched, `ok` | unchanged |
| 3 | `consistent_with` with refs, ceiling `associated` | `overclaimed`, `allowed_claim_level = associated` | recorded as `associated`, reported as a downgrade |
| 4 | `causal` with valid refs | `overclaimed` | still `overclaimed` (never downgraded into silence) |
| 5 | `consistent_with` with no refs | `overclaimed`, `allowed_claim_level = observed` | recorded as `hypothesis`, reported as a downgrade |
| 6 | unsupported ceiling | `overclaimed` | recorded as `hypothesis` |

The general requirement: **`claim_level` on a recorded finding is a function of
the evidence, not of the model's word.** Concretely, the same model result
recorded against two objects with different evidence must yield different claim
levels, and recording it against an object with no evidence must yield
`hypothesis`.

Plus the standing gates: `make check` is 0 errors / 0 warnings / 0 notes, and
the 18 existing tests in `test-ai-claims.R` still pass or are updated with the
reason for each change recorded.

## 4a. Correction found during implementation: the derivation was unreachable

The first implementation was faithful to section 2.2 as written, and is
well-built: it downgrades correctly, keeps `estimated`, and preserves the
`causal` ban. But it **cannot upgrade**, and therefore **has no effect on any
real run**. Measured across every claim level on an object whose ceiling is
`associated`:

| requested | derived | from_ceiling | status |
|---|---|---|---|
| `observed` | `observed` | TRUE | ok |
| `measured` | `measured` | TRUE | ok |
| `associated` | `associated` | TRUE | ok |
| `consistent_with` | `associated` | TRUE | **downgraded** |
| `suggestive` | `associated` | TRUE | **downgraded** |
| `estimated` | `estimated` | FALSE | ok |
| `hypothesis` | `hypothesis` | **FALSE** | ok |
| `causal` | `hypothesis` | FALSE | **overclaimed** |

The derivation runs only inside the `sclet_ai_findings_needing_support` branch
(`consistent_with`, `suggestive`, `causal`). `hypothesis` is not in that set, and
the "no refs" early-return fires before the ceiling is consulted, so a
`hypothesis` finding that cites valid evidence never asks what that evidence
supports.

Since the model emits `hypothesis` for **114 of 114 findings measured**, the
change is a no-op on real data -- section 1's problem, that `claim_level` is a
constant, is not solved.

This is a defect in this spec, not in the implementation. Section 2.2's step 3
said "otherwise", but in the pipeline order the "otherwise" was unreachable for
the one level the model actually uses. The requirement the spec should have
stated plainly:

> **A finding that cites evidence whose ceiling exceeds `hypothesis` must be
> recorded at the ceiling, even when the model asked for `hypothesis`.**
> `hypothesis` is the floor for an *unevidenced* finding, not a cap for an
> evidenced one.

The upgrade is safe for the reason verified in section 5: the ceiling is the
minimum over the supporting nodes' own declared `claim_level`s, so a derived
level can never exceed what the cited evidence already claims for itself.

Amended rule for step 3:

1. `estimated` -> keep (provenance, not strength).
2. `causal` -> problem, `refused = TRUE`, derived `hypothesis`. Never silent.
3. **No refs** -> `hypothesis`, `derived_from_ceiling = FALSE`.
4. **Refs present** -> consult the ceiling for *every* requested level,
   including `hypothesis`, and derive
   `min(requested mapped through sclet_ai_finding_allowed_level(), ceiling)`,
   where an unmapped requested level ranks at 0 (`hypothesis`) so the ceiling
   wins. This is what makes the derived level a function of the evidence rather
   than a copy of the model's word.
5. `overclaimed` is then reserved for the `causal` ban and for unresolvable
   refs -- not for a mere downgrade, which is a correction.

## 4b. Consequence: `overclaimed` narrowed, and what that changes

With derivation in place, `overclaimed` is reserved for the two cases that are
genuine violations -- the `causal` ban, and unresolvable evidence refs. A mere
downgrade is a *correction*, so it no longer stops `RecordAIResult()`.

That changes one existing behaviour, deliberately and with the consequence
verified rather than assumed:

| finding | before | after |
|---|---|---|
| `consistent_with`, no refs | **refused** (`overclaimed`) | recorded at `hypothesis`, `requested_claim_level = consistent_with`, `ceiling_reasons = "no evidence_refs supplied"` |
| `consistent_with`, refs below ceiling | **refused** | recorded at the ceiling, `requested_claim_level` preserved |
| `causal` | refused | **still refused** |
| unresolvable evidence ref | refused | **still refused** |

The invariant is preserved: nothing over-asserted reaches the ledger. An
unevidenced strong claim now lands at `hypothesis` -- the level that asserts
nothing -- instead of blocking the whole result from being recorded. That is a
better outcome than the old one: a refusal discarded the entire result, so the
over-claim was lost rather than corrected, whereas the recorded finding now
carries both the corrected level and the model's original request.

The one thing to keep in view: `RecordAIResult()` no longer refuses an
*unevidenced* claim at all, so `audit_claims = TRUE` is no longer a strict gate
against recording one. It is a gate against recording one *at strength*. If a
future consumer treats "recorded" as "vetted", that reading is now wrong; the
recorded level is the vetted value, and `requested_claim_level` shows the gap.

## 4c. Measured outcome: necessary, not sufficient

Four real runs on pbmc3k after this change landed (`sclet` `a0a4c82`; records in
`sclet-benchmark-ai/results/derived_claim_levels.jsonl`) produced **zero
upgrades**, including in the two conditions built to exercise the upgrade path.
The `with_evidence` condition carried `ev:findallmarkers` with the ceiling
`allowed` at `associated`, verified in the record, so a finding citing it would
have been raised. None did.

Checked directly against a live call: on an object carrying
`ev:findallmarkers`, the model emitted **0 findings with `evidence_refs`** out of
12 and 19, and every finding at `hypothesis`.

So the diagnosis in section 1 was incomplete. It said the ladder was inert
because the model always takes the default level. The model also supplies no
evidence references -- and the ceiling reasons about *evidence*, so without a
reference there is nothing to reason about and the derivation correctly lands
every finding at `hypothesis` through its "no refs" path.

This spec therefore fixes the right half and cannot, on its own, make
`claim_level` informative. That is a property of the result, not a defect in the
implementation, and it is worth stating plainly rather than leaving the reader
to infer that the change underperformed.

The remaining question is open, and section 3's rejection still applies: do not
instruct the model to cite evidence and claim more strongly, because its current
non-fabrication is a measured virtue and that instruction trades it for a
number. The plausible routes are (a) have the engine attach the evidence a
finding relates to, since it knows which evidence each step created, leaving the
model only to describe; or (b) accept that `claim_level` stays at the floor for
interpretive reviews and stop treating the field as informative. Both need their
own spec.

**Closed for route (a), and now measured online.**
`.dev/spec-ai-evidence-authorship.md` implements it: the engine attaches the
evidence a finding relates to when the mapping is forced (sole completed
non-decision node), and declines otherwise. Offline: a `hypothesis` finding
that cited nothing on a one-node object is recorded at that node's own level
(`associated`) with `evidence_refs_attached = TRUE` and the rule name recorded
beside it; on a two-node object the floor holds and the refusal is reported.
Online (2026-10-13, `sclet-benchmark-ai` `derived_claim_levels.jsonl`): the
zero-upgrade rows above became 6/7, 12/13 and 9/10 upgrades on the
`with_evidence` condition, while the `without_evidence` control still lands
every over-requested finding on `hypothesis`. `claim_level` is no longer a
constant -- it varies with what the ENGINE can attach, not with what the model
claims, which is exactly the asymmetry §3's rejection demanded. The remaining
question the online data raises is the one this spec's rule settled offline
and the real model re-opened: requests BELOW the ceiling are raised to it, so
a finding that self-reports `observed` for a structural fact gets recorded at
`associated`. See section 9 of the authorship spec.

## 5. Risk

The main risk is **over-correcting in the other direction**: if the engine
derives `associated` from a single `ev:findallmarkers` node, recorded findings
will get a stronger-looking level than the model chose, and a reader may take
that as the system endorsing a claim. Two mitigations, both required:

- always store `requested_claim_level` alongside `claim_level`;
- have `sclet_ai_claim_ceiling()`'s `reasons` (already returned: "all supporting
  evidence resolves to 1 independent line(s) of support; a claim cannot exceed
  'associated'") travel with the recorded finding, so the level is never quoted
  without its basis.

A second, smaller risk: `n_independent_groups = 1` means the ceiling is
`associated` for essentially every single-node object, so the derived level will
still be near-constant in practice -- just at a level that carries information
rather than at the floor. That is an improvement, not a solution to the
"one node is thin evidence" problem, and should not be over-sold.

### Why the upgrade path is bounded anyway (verified)

Deriving can move a finding *up* from `hypothesis` to `associated`, which is the
change most likely to be misread as the system endorsing a claim. It is bounded
by construction, and this was verified in the implementation rather than
inferred from intent:

- the ceiling is a minimum over the supporting nodes' own `claim_level`s -- the
  implementation states this in its own reason string, "bounded by the weakest
  supporting node" (`R/ai-claims.R`);
- measured on a pbmc3k object carrying `ev:findallmarkers`: that node declares
  `claim_level = associated`, and the ceiling independently computes to
  `associated` with `n_independent_groups = 1`.

So a derived level can never exceed the strength the cited evidence already
declares for itself. An upgrade to `associated` means "this finding is no
stronger than the evidence it cites" -- not "the system endorses a biological
conclusion". Wherever a derived level is surfaced, that reading should be stated,
because the bare number does not convey it.
