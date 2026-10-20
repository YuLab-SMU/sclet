---
name: sclet-claims
description: >
  Claim-calibration guidelines for sclet's interpretive entry points
  (AIReviewQC, AIExplainAnalysis, AIStatus and friends): the claim-level
  ladder and what each level is allowed to assert, how evidence_refs and the
  claim ceiling bound what can be written, the calibrated language for each
  level, and the fabrication traps specific to single-cell data (naming
  genes that are not on the object, reading clusters as cell types,
  treating a null result as a negative finding). Use when writing findings,
  answering interpretive questions about an object, grading or auditing
  claims, or reporting QC, marker, annotation or trajectory results.
---

# Claim-calibration guidelines for sclet findings

These guidelines implement the calibration tables proposed for single-cell
transcriptomics interpretation ("From Cluster to Claim: Calibrating
Interpretation in Single-Cell Transcriptomics", DOI 10.1002/ggn2.70045) on
top of sclet's own claim machinery (R/ai-claims.R: the finding claim levels,
the evidence ceiling, and the deterministic claim auditor). The engine
enforces the ceiling; these guidelines are about producing claims that are
honest BEFORE the engine grades them, and about the failure modes specific
to single-cell data.

## The ladder, and what each level may say

| level | asserts | allowed language | requires |
|---|---|---|---|
| `hypothesis` | nothing; explicitly speculative | "could", "would be consistent with", "speculatively" | nothing |
| `observed`/`measured` | a structural fact read off the object | "the object has", "is recorded as" | the fact is on the object |
| `associated` | a statistical association between variables | "is associated with", "differs between" | >= 1 independent evidence line |
| `consistent_with` | evidence supports a specific biological reading | "is consistent with", "suggests" | >= 2 independent evidence lines |
| `causal` | mechanism | — | **never produced by the AI** |

`causal` is banned outright: a single-cell analysis does not establish
mechanism, and the system refuses rather than downgrades it. If an
explanation feels causal, write it as `consistent_with` and name the
alternative explanations.

## The ceiling is a function of the evidence, not of confidence

Recorded levels derive from cited evidence, bounded by the weakest
supporting node:

- a finding citing no evidence lands at `hypothesis` no matter how it is
  phrased;
- one independent line of support caps at `associated`;
- two or more independent lines can reach `consistent_with`;
- re-describing one measurement twice (a marker test and the enrichment
  derived from it) is ONE line of support -- the ledger tracks parentage so
  corroboration is not invented.

Independent means independent: a second engine over the same expression
matrix and the same group labels is not a second measurement. This is the
corroboration trap the dependency groups exist for.

## Single-cell fabrication traps

These are the specific ways an interpretive answer over-reaches on
single-cell data. Each has been measured in this system's benchmark
corpus:

1. **Naming genes that are not on the object.** Gene names are the most
   seductive hallucination in this domain because the vocabulary is dense
   and plausible (`CD8A`, `MS4A1`). Only name a gene that is in the object's
   rownames; when citing a marker, cite the evidence node that measured it.
2. **Reading clusters as cell types.** A cluster is a group of cells under
   one resolution at one parameter set. Without annotation evidence it is
   `cluster_N`, not `T cells`. An annotation transferred from a reference
   carries that reference's evidence; a cluster marked by two markers is
   provisional.
3. **Treating a null result as a negative finding.** "No significant markers"
   on an object with two near-identical groups is the correct result, not a
   failure to find them. Report what the data supports -- including that the
   groups are not separable -- rather than promoting borderline genes to
   fill the expected shape of an answer.
4. **Reading a proportion shift as expression change.** A cluster shrinking
   between conditions is a composition change; without a differential-
   abundance test it is not DE, and implying it is, misstates the biology.
5. **Promoting a QC observation to a population claim.** High mitochondrial
   fraction in one group is a QC fact; "dying cells form their own cluster"
   is a biological claim that needs marker evidence.
6. **Inventing proposed actions.** `proposed_actions` carries registered
   action names; a prose sentence ("compute QC metrics and inspect them")
   belongs in `recommendations`, not there.

## Structural facts are provenance, not findings

"The object has 11,769 cells", "clustering ran at resolution 0.8 with seed 1",
"the KEGG run used ENTREZID" -- these are facts read off the object and its
records, and they belong in the report's methods/provenance section, not in
`findings`. A finding with no `evidence_refs` derives to `hypothesis` (the
engine cannot verify a statement it has nothing to check), so writing a
structural fact as a finding either downgrades a true statement to a
hypothesis or invites inventing a citation for it. The renderer prints the
recorded analyses and the executed steps for exactly this purpose; cite a
record only when the claim is *supported by* it.

## Honest reporting of refusals and missing evidence

- If the evidence ceiling is `absent` (no evidence nodes), say what the
  object does NOT carry and what would create it (executing a plan with a
  named step). Do not treat it as permission to speculate freely.
- If a gate refused an action (no annotation reference, no design
  confirmation, missing modality), report the refusal and the unblocking
  step. Never report the substituted analysis as the answer.
- If findings are trimmed or levels downgraded by the auditor, the record
  keeps `requested_claim_level` -- do not re-assert the requested level in
  prose after the engine lowered it.

## The floor is a virtue

Measured across this system's real-provider corpus: on a null object the
model named zero genes and asserted zero claims above the floor, and every
finding sat at `hypothesis` unless evidence supported more. That is the
correct behaviour and it is worth more than a confident-sounding answer:
on single-cell data the cost of an over-claim is a fabricated biological
narrative, and the cost of `hypothesis` is one clarifying sentence.
