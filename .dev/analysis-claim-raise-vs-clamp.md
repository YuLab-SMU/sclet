# Decision memo: should the derived claim level ever RAISE a request?

Date: 2026-10-13. Read-only analysis of `R/ai-claims.R` @ `devel`; no code changed.
**Recommendation: (c) clamp-to-ceiling + raise-from-floor, PLUS the relevance guard
`spec-ai-evidence-authorship.md` §4 req. 2 already requires and the code never
implemented.** One line: the raise is the only thing authorship exists for and the clamp is
the only thing the ladder ordering demands — but the Lead's C5 is not on the raise/clamp
axis, and the only rule on that axis that stops it is (b), which stops it by disabling the
mechanism.

## 1. The mechanism, with the ladder arithmetic

Evidence ladder (line 1): `observed(1) < measured(2) < associated(3) <
consistent_with(4)`; `hypothesis` and `estimated` are absent, so `sclet_ai_claim_rank()`
returns `NA` for both (verified). §3c (R/ai-claims.R:283-298) ignores the request's rank
entirely: `derived <- if (ceiling_rank < 1L) "hypothesis" else
sclet_ai_claim_from_rank(ceiling_rank)` — a pure function of the ceiling, so any request
strictly below it is overridden upward. Path: the DE node is hard-coded `claim_level =
"associated"` (R/ai-execution.R:939); on a one-node object the ceiling is
`min(3L, min_rank)` = `associated` with `n_independent_groups = 1`, status `allowed`; the
model's `observed` (rank 1) becomes `associated` — "I looked, this is what is there"
rewritten into a claim that two variables co-vary.

Upgrades verbatim, `derived_claim_levels.jsonl` last 6 rows (2026-10-13, `deepseek-chat`,
pbmc3k): with_evidence/qc_review `observed -> associated` x6 of 7; with_evidence/
explain_analysis `observed -> associated` x9 + `measured -> associated` x3 +
`hypothesis -> associated` x1 of 13; with_evidence/status `observed -> associated` x9 of
10. without_evidence (3 tasks): 0 upgrades, `observed -> hypothesis` x20 and
`measured -> hypothesis` x10. `measured -> associated` is the same defect one rank lower,
so the harm is broader than "a different kind of claim".

## 2. Candidate rules on six cases

C5 is the Lead's case, reproduced offline (§7): sole `associated` node, cite-nothing
`hypothesis` finding whose text is "no mitochondrial or ribosomal content check and no
doublet detection were recorded" — a claim that evidence is ABSENT. Rule
`sole_evidence_node`, recorded `associated`.

| rule | C1 `observed`, ceiling assoc, cited | C2 `hypothesis` + attached node | C3 `consistent_with`, ceiling assoc | C4 no level, no evidence | C5 `hypothesis` + attached IRRELEVANT node | C6 `hypothesis` + MODEL-cited node |
|---|---|---|---|---|---|---|
| (a) ceiling-only (today) | **associated** | associated | associated | hypothesis | **associated** | associated |
| (b) `min(requested, ceiling)` | observed | **hypothesis** | associated | hypothesis | hypothesis | hypothesis |
| (c) clamp; raise from floor | observed | associated | associated | hypothesis | **associated** | associated |
| (d) clamp; raise if ENGINE attached | observed | associated | associated | hypothesis | **associated** | hypothesis |
| (e) clamp; raise if MODEL cited | observed | hypothesis | associated | hypothesis | hypothesis | associated |

**C2 and C5 are the same row for the derivation.** Its only inputs are `object` and the
finding; post-attachment both are `requested = hypothesis` with one `associated` ref.
Nothing in its input separates the legitimate raise from the false one.

To the Lead's question, directly: **(b) does stop C5**, via the same arithmetic that stops
C2 — `min(0, 3) = 0`. That is not a coincidence to bank on; it is blanket suppression of
every floor request, which is authorship's entire payload (§4c: 0 of 23 and 0 of 13 moved
before it). On a one-node object — the only object authorship fires on — (b) makes
authorship a pure no-op. (d) is the worst option: keeps C5, drops the vetted case. **(e) is
the only rule on the axis that separates C2 from C5 without a text check**, because
`attachment_rule` is computed one line earlier in the same loop (R/ai-claims.R:360-362) and
could be passed down.

## 3. The real C5 defect, and why (c) + guard

The defect is upstream of the derivation. Authorship spec §4 req. 2 permits attaching a
sole node only if it "names a step the finding's text does not contradict". **That clause is
not implemented:** `sclet_ai_attach_finding_evidence()` (lines 97-153) filters on completed,
non-`user_decision`, exactly-one, and never reads the finding's text. The Lead's numbers on
`interpret_evidence_refs.jsonl` records 3-4 confirm the exposure is the attached population,
not the cited one: of **51** cite-nothing findings, **46** never name the node anywhere in
their prose (verified by splitting `findings_flat`), and every "no X was recorded" finding is
in the 46. Attachment with no relevance test is the bug; the raise only propagates it.

- **(c)** answers the question actually posed and is the minimum change killing
  `observed -> associated` and `measured -> associated` while keeping the raise.
- **(b) rejected**: does not fix relevance, disables the mechanism.
- **(e) is the interim if the guard is refused** — the only option that stops C5 and keeps
  the model-cited raise, needing no attribution and no text heuristic. Its cost: it reverts
  the three assertions that ARE authorship's headline (§7), i.e. it decides "attachment
  records provenance, attachment does not confer strength" — defensible, but a reversal of
  that spec's §1, not an amendment.
- **Strongest argument against (c) + guard**: the guard is the hard part and (c) does only
  the easy part. A text-contradiction check is either an unpredictable heuristic or a model
  call, and §3 of the derived-levels spec rejects leaning on the model. If the guard is never
  written, (c) ships C5 unchanged and (e) would have been the better choice. Stating that
  rather than assuming the guard lands.
- **Second, independent objection to (c)**: it makes the recorded level partly a function of
  the model's self-report again — what §2 of the derived-levels spec set out to end — and
  that self-report is uncalibrated: the same batch volunteers `measured` x10 on an object
  with **zero** evidence nodes, where `measured` is false. The no-refs floor saves those, but
  on a multi-node object (attachment declines, clamping fires) a `measured` guess is recorded
  at `measured`. (c) trusts direction while the ceiling controls magnitude; (b) and (e) trust
  it more, not less. Only (a) refuses the self-report, and (a) manufactures over-claims.

## 4. Does clamp-only (b) resurrect the 114-of-114 problem? Yes.

`hypothesis` must rank 0 (raw `min()` on `NA` returns `NA`; `sclet_ai_claim_from_rank(NA)`
returns `NULL`, which would write an empty `claim_level`), so `min(0, ceiling_rank) = 0` →
`hypothesis`. Authorship still attaches, the ceiling still computes `associated`, the record
still reads `hypothesis`: the attached ref buys nothing.

Does attachment make `min()` safe now when it was not then? Partly — the distinction is
which branch fires. At the 114/114 measurement (2026-10-12, first 4 jsonl rows,
`hypothesis=23/13/17/17`) authorship did not exist, so every finding took the **no-refs**
branch and `min()` was never reached; pre-§4a had a second, earlier inertia (derivation ran
only inside `sclet_ai_findings_needing_support`, which excludes `hypothesis`). Today a
cited-nothing `hypothesis` finding reaches the ceiling because the engine gave it a ref — but
reaching the ceiling is not being raised by it. So bare `min()` still cannot be adopted;
`min()` with a floor exception can, and that is (c). Not "attachment exists, therefore
`min()` is safe", but "attachment exists, therefore the floor must be exempted from `min()`
or the attachment is decoration". (e) conditions the exemption on who supplied the ref.

## 5. Vocabulary asymmetry

- **`estimated`** returns at lines 185-194 *before* refs are read (`ceiling_status =
  "provenance"`) and attachment declines it (line 113), so no ref is ever attached. It is
  rank `NA` — the trap: if the early return moved after the rank comparison, `estimated`
  would rank 0, look like the floor, and be raised to `associated`, rewriting provenance as
  strength. **No candidate rule harms it provided that guard stays step 1.**
- **`suggestive`** maps down to `consistent_with` (line 66) before the ceiling, so under
  (c)/(e) it ranks 4 and is only clamped, never raised.
- **`causal`** — banned, `refused = TRUE`, never a silent downgrade; orthogonal to all five.
- Dead code: `sclet_ai_findings_needing_support` (line 58) is referenced nowhere in `R/` or
  `tests/` (grep-verified); a leftover of the pre-§4a gating, not live behaviour.

## 6. Attribution: the Lead's falsification holds, and it removes my first instinct

Verified offline (§7): `sclet_ai_target_matches("ai_evidence", "ev:findallmarkers", t)` is
`TRUE` for `"markers"`, `FALSE` for `"run_de_test"`, `"findallmarkers"` and an execution id.
The matcher greps `tolower(paste(type, id))` (R/ai-context.R:126-154) against a fixed alias
table; the node's `source` (R/ai-execution.R:942) is never part of the matched string. So
**scoping attachment by `target` is falsified and I do not lean on it**: (c) needs the text
guard, (e) needs `attachment_rule`, computed locally in-process — neither needs attribution.
Any future proposal to promote attached findings by execution provenance must first tag
evidence nodes with the producing execution id and change what the matcher greps: a
data-model change, not a filter.

## 7. What would have to change (nothing changed here)

Code:
- `sclet_ai_derive_finding_claim_level()` §3c: `mapped <-
  sclet_ai_finding_allowed_level(requested)`; `req_rank <- sclet_ai_claim_rank(mapped)`,
  `NA -> 0L`; `final_rank <- if (req_rank < 1L) ceiling_rank else min(req_rank,
  ceiling_rank)`. For (e), take an `attachment_rule` argument and use
  `if (req_rank < 1L && identical(attachment_rule, "model_cited")) ceiling_rank else min(...)`.
  Rewrite the §3c comment: "deliberately NOT `min(requested, ceiling)`" (lines 289-295)
  stops being true, and the 114/114 justification retires with a note that the model now
  self-reports `observed`/`measured`.
- `sclet_ai_attach_finding_evidence()`: implement spec req. 2's text guard — the actual C5
  fix; at minimum decline attachment when the finding asserts absence.
- `sclet_ai_audit_result_claims()`: no logic change, but the label is wrong today — an
  UPGRADE gets `entry_status = "downgraded"` and increments `n_downgraded` (verified:
  `observed -> associated` → entry "downgraded", `n_downgraded = 1`, overall "downgraded").
  That is why the jsonl reads `audit_status = "downgraded"` on rows with 12 upgrades and 0
  downgrades. Split `raised`/`corrected` or (c) makes the field read `ok` on exactly the
  rows that should show a correction.
- Unchanged: `sclet_ai_claim_ceiling()`, `sclet_ai_recorded_findings()`, `RecordAIResult()`.

Tests. **Nothing today covers C1** — which is why the over-claim shipped;
`test-ai-evidence-authorship.R:292-298` tests a low *ceiling*, not a low *request*. Under
(c), `test-ai-derived-claim-levels.R:204` ("the request is raised to the ceiling", asserting
`consistent_with -> associated`) and `test-ai-claims.R:100-133`, `:156-174` all stay green —
the first is mis-*titled*, not wrong. Under **(e) three load-bearing authorship assertions
fail**: `test-ai-claims.R:100-133` (line 122), `test-ai-claims.R:156-174`,
`test-ai-derived-claim-levels.R:88-110` (line 105), plus that spec's §5 criterion 1. Either
rule needs new cases C1-C6 asserted; (c) also needs `measured -> associated` asserted as NOT
an upgrade.

Specs: `spec-ai-derived-claim-levels.md` §4a step 4 prescribes `min(requested mapped,
ceiling)` "where an unmapped requested level ranks at 0 ... so the ceiling wins" —
self-contradictory as written (`min(0,3) = 0` is the floor). The code resolved that toward
always-ceiling; **the spec text is closer to (c) than the code is**, so (c) is the faithful
reading and (a) the over-correction. §5's risk list already named this failure mode.
`spec-ai-evidence-authorship.md` §4 req. 6 pins "still not `min(requested, ceiling)`" and
req. 2's guard is unimplemented; both must be reopened by this memo.

## 8. Verification (offline, no provider, no billed call)

```bash
cd /home/wang/data/source/omics/sclet && Rscript -e '
suppressMessages(pkgload::load_all(".", export_all=TRUE, helpers=FALSE, attach_testthat=FALSE, quiet=TRUE))
sce0 <- SingleCellExperiment::SingleCellExperiment(list(counts=matrix(1,3L,4L)))
node <- function(o, id, claim="associated", src="s", dg="g1", kind="deterministic_summary")
  RecordAIEvidence(o, list(id=id, kind=kind, values=list(f=1L), claim_level=claim, source=src,
  dependency_group=dg, parents=character()), source=NULL, parents=character(), scope=NULL)
sce1 <- node(sce0, "ev:findallmarkers", src="findallmarkers", dg="g_fm", kind="test")
a1 <- function(sce, f) sclet_ai_audit_result_claims(sce, list(findings=list(f)))$findings[[1L]]
show <- function(l, e) cat(sprintf("%-32s -> %-11s req=%-15s rule=%-18s status=%s\n", l, e$claim_level,
  e$requested_claim_level, e$attachment_rule, e$status))
show("C1 observed + ceiling assoc", a1(sce1, list(claim_level="observed", evidence_refs="ev:findallmarkers")))
show("C2 hypothesis + attached node", a1(sce1, list(claim_level="hypothesis")))
show("C3 consistent_with + ceiling", a1(sce1, list(claim_level="consistent_with", evidence_refs="ev:findallmarkers")))
show("C4 no level + no evidence", a1(sce0, list()))
show("C5 absent-evidence claim", a1(sce1, list(claim_level="hypothesis",
  statement="no mitochondrial or ribosomal content check and no doublet detection were recorded")))
show("C6 hypothesis + model-cited", a1(sce1, list(claim_level="hypothesis", evidence_refs="ev:findallmarkers")))
cat("ranks:", sclet_ai_claim_rank(c("observed","measured","associated","consistent_with")),
  "| hypothesis NA:", is.na(sclet_ai_claim_rank("hypothesis")),
  "| from_rank(NA) NULL:", is.null(sclet_ai_claim_from_rank(NA)), "\n")
print(sclet_ai_claim_ceiling(sce1, evidence_ids="ev:findallmarkers")[c("ceiling_claim_level","status","n_independent_groups")])'
```

Expected today (a): `associated, associated, associated, hypothesis, associated,
associated`. Under (c): C1 → `observed`, C5 still `associated` (unfixed — hence the guard).
Under (e): C1 → `observed`, C2/C5 → `hypothesis`, C6 → `associated`.

Attribution falsification (§6):

```bash
cd /home/wang/data/source/omics/sclet && Rscript -e '
suppressMessages(pkgload::load_all(".", export_all=TRUE, helpers=FALSE, attach_testthat=FALSE, quiet=TRUE))
for (t in c("markers","run_de_test","findallmarkers","ai_exec_20261013")) cat(sprintf("%-18s -> %s\n", t,
  sclet_ai_target_matches("ai_evidence", "ev:findallmarkers", sclet_ai_normalize_target(t))))'
```

51/46 counts (§3): records 3-4 of `interpret_evidence_refs.jsonl`, split `findings_flat` on
`" | "`, `evidence_refs=` empty → 51; of those, no `findallmarkers` substring in the text →
46. Dead helper: `grep -rn sclet_ai_findings_needing_support R/ tests/` → `R/ai-claims.R:58`
only.
