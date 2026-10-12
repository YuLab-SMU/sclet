# Spec: derive finding claim levels from evidence instead of asking the model

Status: proposal (not implemented)
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
  level into the recorded finding, never the requested one. That is what makes
  the ledger's claim levels meaningful for anything reading them later.
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
