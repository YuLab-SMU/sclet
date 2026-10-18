# Malignant-cell identification: from the iterative port to a two-line call

Status: implemented in `feat/malignant-cnv` (worktree `../sclet-malignant-cnv`).
Upstream discussion: [YuLab-SMU/ProjectYulab#26](https://github.com/YuLab-SMU/ProjectYulab/issues/26).

## What the issue is actually about

The issue asks for an implementation of the malignant-epithelium classification
from Zhang et al., *Gut* 2020 (doi:10.1136/gutjnl-2019-320368): bulk tumour
versus normal DEGs -> score single cells -> k-means into two -> re-derive DEGs
from the split -> repeat until the labels stop moving.

The third comment is the valuable part. It reproduces the method and then argues,
from the paper's own supplementary tables, that the "malignant" signature is
mostly intestinal-function genes (CLDN4, CLDN7, TFF3, REG4) and that what the
iteration separates is **gastric identity versus intestinal identity**, not
malignant versus benign. Two further points:

* **Sensitivity.** Changing only the seed gene source moved the malignant rate
  from 30.6% to 6.6%.
* **Reference drift.** After the first round the "non-malignant" pool contains
  tumour-sample cells, so the comparison baseline is no longer the normal
  samples. The paper's own results describe SPEM in the "no malignant" group.

All three of these are structural consequences of the algorithm, not
implementation bugs:

1. K-means on two scores always returns two groups, including when no two-way
   structure exists.
2. Re-deriving the signature from the split is a positive feedback loop: each
   round amplifies the axis that separated the last round. The axis that
   maximises variance need not be the malignant axis.
3. The reference is endogenous. It is whatever the split produced, so it moves
   with the split.

## What we built

### Layer 1: the reference port (`R/malignancy.R`)

`FindMalignantSeedGenes()` and `ClassifyMalignantEpithelium()` implement the
paper's procedure faithfully, as a baseline to compare against. Three failure
modes were found and handled while writing it:

| Failure | Consequence if unhandled | Handling |
|---|---|---|
| `AddModuleScore` control bins | Scoring a signature-only submatrix makes signature genes each other's controls; every score collapses to 0 and k-means has nothing to split | Score against the full feature set; exclude both signatures from each other's control pool |
| k-means always splits | A 55/5 fixture came back 92/8 and would have been reported as a classification | `split_degenerate` stop reason plus a warning; never silently returned |
| Flipped cluster ids | The raw diagonal reads a mirrored assignment as 0% agreement, so the loop never converges | Agreement = `max(same, swapped) / n` |

The port reports diagnostics rather than a verdict: iteration log, side-of-origin
composition (showing the reference drift explicitly), seed-gene tracking for a
competing identity programme, and a standing caveat that converged labels are
consistent, not validated.

### Layer 2: CNV evidence (`R/malignancy-cnv.R`)

Malignancy's defining feature is genomic instability, which is orthogonal to any
transcriptional identity programme. `RunCNVInference()` runs inferCNV against a
**caller-supplied, exogenous reference** -- normal cells or normal samples named
by the user, never cells the algorithm itself classified. The reference therefore
cannot drift.

`check_cnv_readiness()` names every missing prerequisite. When inferCNV or its
JAGS system library is absent, `RunCNVInference()` records
`status = "unavailable"` with the reason and writes `NA`, **not zero**: an
absent measurement must never read as a negative result.

Implementation note that cost real debugging time: inferCNV's `expr.data` holds
reference-centred smoothed values whose level sits near 1 for every cell.
`colMeans(abs(expr))` therefore measures the baseline, and on inferCNV's own
bundled oligodendroglioma example the malignant and normal groups differ by
**-0.005**. Subtracting each cell's own median first leaves the profile shape,
which separates the same groups by **+0.038** (malignant 0.088-0.108, reference
0.051-0.058). A median, not the reference mean, so that a large amplified tract
cannot pull its own baseline toward itself.

### Layer 3: the generic evidence registry (`R/phenotype-evidence.R`)

Extracted from the malignant call rather than bolted on afterwards: the first
version hard-coded inferCNV as "the second line", which was a design defect, not
just a missing generalisation.

`PhenotypeEvidenceProviders()` holds named providers. A provider is a plain list
(`available()`, `decide()`, `dependence`), so a user registers their own without
subclassing anything. `CallPhenotypeCells()` counts lines and reports a four-state
call: `supported` (two or more positive), `probable` (one), `undetermined`
(conflict, nothing usable, or a competing programme explains the split),
`refuted` (an explicit negative and no positive).

Three rules keep it honest:

* **`dependence`.** A line computed from the same expression matrix as the
  signature is `shared_with_signature` and cannot count toward `supported`. This
  is the property the iterative method in `R/malignancy.R` lacks.
* **Missing is not negative.** A low signature score means the line did not fire,
  not that the phenotype is absent; only a provider that says `negative` counts
  as a no. A signature-negative cell with one positive provider is `probable`,
  not `undetermined` -- and `refuted` requires an explicit negative.
* **Claims reuse the package's ladder** (`observed` -> `measured` ->
  `associated` -> `consistent_with`) instead of inventing a second vocabulary for
  confidence. Nothing reaches `observed`: a per-cell label is never a direct
  observation of cell state.

The iterative search is deliberately absent here. Iteration is appropriate when
the cell population itself is unknown; a call built on named evidence is not
improved by re-deriving its own signature.

### Layer 4: the four-state malignant call (`R/malignancy-combine.R`)

`CallMalignantCells()` is now thin: it names the CNV provider, maps the generic
call onto the malignant vocabulary, and requires a CNV record rather than
silently degrading to a signature-only call. It refuses to resolve disagreeing
evidence by majority.

| State | Condition |
|---|---|
| `malignant` | CNV-positive **and** signature-positive |
| `likely_malignant` | exactly one line positive |
| `undetermined` | conflict, CNV unavailable, or a competing identity programme explains the split |
| `non_malignant` | both negative |

Two properties are pinned by tests:

* With CNV unavailable, **no cell reaches `malignant`**. Signature-positive cells
  become `likely_malignant`; the rest `undetermined`.
* A CNV cut derived from the reference population must sit strictly above that
  population's maximum; otherwise a constant or exactly-measured reference sits
  on its own threshold and reads as CNV-positive (mad = 0).

The `identity_genes` argument makes the issue's central criticism operational:
score a suspected competing programme (for gastric cancer, the intestinal-
metaplasia profile), and if it separates the two evidence-positive sides beyond
`identity_demote`, demote those cells to `undetermined` and say so.


## What the audit found (and why these were not test failures)

The first round of tests all passed while two of the framework's central promises
were broken. Neither defect was a boundary case: one disabled the competing-
hypothesis test, the other disabled the independence rule.

**Circular competing hypothesis.** The test used `call` to define "the cells
under discussion". Since the signature selects those cells, the competing
programme was guaranteed to separate them -- the question answered itself. Worse,
the demotion swept `refuted` cells into `undetermined`, so a refusal could be
overturned by the very test meant to challenge the affirmative calls. Fixed:
groups are now defined by the evidence (`in_scope` = at least one line positive,
`residual` = none), never by the call, and `refuted` is never demoted.

**Independence by declaration only.** Registering one measurement under two names
made every cell `supported`. The registry cannot see inside a provider, but it
can honour a declaration: providers sharing a `source` now count as one line, and
re-registering an existing name is refused instead of silently replacing it.

Three smaller ones, also fixed: `n_cnv_regions` was written to colData and the
state record but never computed (always `NA`, a promised measurement that did not
exist); a non-count assay was coerced with `mode(x) <- "integer"`, truncating
every value silently; and the generic-to-malignant mapping had two overlapping
rules that agreed only by evaluation order.

The lesson worth keeping: the defects were all in the *meaning* of the columns,
and "every assertion passes" said nothing about them. Each fix now has a
regression test that fails without it.


## What the second audit found

A second pass, reading the implementation rather than the tests, found two
crash-level defects and three silent ones. The pattern repeats: **every defect
was in how a value is obtained, not in how it is used.**

**A gene order with no overlap crashed.** Normalising the `chr` column happened
*after* filtering to matching genes, so an order sharing no id with the object
assigned a column to a zero-row data frame: `replacement has 1 row, data has 0`.
Gene-id mismatch (ENSG vs symbol) is the ordinary way to hit it. Normalisation
now precedes subsetting, and no overlap is a named error that prints both id
styles.

**The readiness gate counted genes but never checked chromosomes.** inferCNV
slides a window *along* each chromosome, so 400 genes piled onto chr1 passed a
`nrow >= 100` gate and reported `ready: TRUE`. Chromosome coverage is now its own
check. Duplicated gene rows are merged into one interval per gene and reported as
a note that never blocks the run.

**A provider's cell names were discarded.** A provider returning its states named
by cell, in a different order than colData, was matched positionally: every cell
received another cell's verdict, silently. Named states are now matched by name,
and names that do not match colData are refused.

**Cells dropped by inferCNV vanished.** Its minimum-count filter routinely removes
cells, which then landed in colData as unexplained NA. The count is now reported
in a message, in the analysis record and in the state record.

One ordering lesson worth keeping, because it cost two wrong turns here: an
*input error* and a *missing prerequisite* must not share a gate. A
non-count assay and too few reference cells both used to be reported only after
the readiness check, which returned early for an unrelated reason and swallowed
them. Both are now validated before the gate, so the caller hears about the
matrix before the gene order.

## What this deliberately does not claim

There is no per-cell ground truth for malignant identity. `cnv_cut` and
`transcriptomic_cut` are analysis parameters, not discovered thresholds; moving
them moves the answer. CNV-positivity is strong positive evidence and
CNV-negativity weak negative evidence, because copy-number-quiet but hypermutated
tumours exist and low-depth cells are noisy. Report the state distribution and
the disagreement cells alongside any conclusion.

## Verification

* `tests/testthat/test-malignancy.R` -- 39 assertions, iterative port.
* `tests/testthat/test-malignancy-cnv.R` -- 25 assertions, CNV layer and call.
* `tests/testthat/test-phenotype-evidence.R` -- 30 assertions, registry and
  generic call, including that a shared-matrix provider cannot manufacture a
  two-line call.
* Full suite 653 tests / 0 failures / 5882 passed.
* The burden construction is verified against a real inferCNV object, not only a
  synthetic one; inferCNV itself is mocked where a full run is needed because it
  needs JAGS and takes minutes per small sample.

## Open question for the next round

The issue thread deserves a reply once this is settled (explicitly deferred by
the maintainer to after layers 1 and 2). The reply should state: the port exists
as a baseline, the three failure modes are now named and tested, and the
package's answer separates CNV evidence from transcriptional evidence instead of
re-deriving one from the other.

## First real-data application (2026-10-19)

Data: GSE144735 (KUL3 CRC 10x, raw UMI counts). 300 normal immune/stromal cells
as the exogenous reference, 800 tumour epithelial cells and 300 normal
epithelial cells as observations (1,400 cells, 15,132 genes after filtering).
Transcriptional line derived **per patient** from tumour-vs-normal epithelium
log-fold changes, keeping genes up in at least 80% of patients and dropping
ribosomal/mitochondrial families (IFI6, IFITM3, CEACAM6, BST2, ISG15, CDKN2A,
TFF3, AZGP1, ...). CNV line default: inferCNV, HMM off, reference quantile 0.99.

| measurement | result |
|---|---|
| AUROC CNV burden, tumour vs normal epithelium | 0.709 |
| AUROC CNV burden, tumour vs immune reference | 0.824 |
| AUROC signature score, tumour vs normal epithelium | 0.894 |
| CNV-positive cells at q99 (default) | 5 / 800 tumour, 0 / 300 normal epithelium |
| four-state call at q99 | 1 malignant, 4 likely, 637 undetermined, 158 non-malignant (tumour epithelium) |
| normal epithelium called malignant/likely | 0 / 300 |
| threshold q95 (cut 0.113) | tumour 25.2% malignant/likely, normal epithelium 16.7% |
| threshold q90 (cut 0.102) | tumour 40.8%, normal epithelium 24.3% |

What this says: the signature line works and is biologically plausible, but the
CNV-burden line separates tumour from normal epithelium only moderately, so the
two-line call is specific and mostly **undetermined** rather than overconfident.
Lowering `cnv_cut` buys sensitivity at a steep cost in normal-epithelial false
positives -- the threshold is doing real work, and the default is conservative.
The earlier 300-cell smoke run through the AI action layer (with an illustrative
epithelial marker score) exercised the plumbing only; these numbers are the
first biological read.

Caveat: one dataset, one 1,400-cell subset, in-sample signature, and inferCNV
run without HMM. The next step if sensitivity matters is not to lower the cut but
to improve the CNV evidence (HMM region calls, larger/cleaner reference, or a
different burden statistic) and re-check the same AUROC.
