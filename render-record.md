# Book render record

- Rendered: 2026-10-23 21:11:20 CST
- sclet version: 1.1.0
- sclet library: /home/wang
- sclet commit: b9e996d
- R: R version 4.6.1 (2026-06-24)
- Platform: x86_64-pc-linux-gnu
- SCLET_BOOKDOWN_EVAL: (unset, defaults to true)
- SCLET_BOOKDOWN_EVAL_AI: 1
- aisdk available: TRUE

Counts are per chapter, over three groups:

- **executed** — no `eval=` option, or `eval=TRUE`. Always runs.
- **conditional** — `eval=` is a variable. `eval_chunks` (the
  `SCLET_BOOKDOWN_EVAL` switch) defaults to TRUE, so most of these DO run;
  `eval_ai` defaults to FALSE; `svp_ready` / `pipeline_ready` /
  `file.exists` / `requireNamespace` depend on data this machine may not
  have.
- **off** — `eval=FALSE`. Never runs.

**A record is only comparable to another record taken on a machine with
the same data.** The point of the record is the diff: a chapter whose
numbers drop means something stopped running.

| chapter | chunks | executed | conditional | AI-gated (ran?) | off |
|---|---|---|---|---|---|
| ai-driver.Rmd | 12 | 0 | 12 | 0 (-) | 0 |
| ai-evidence-governance.Rmd | 3 | 0 | 3 | 0 (-) | 0 |
| ai-foundations.Rmd | 3 | 0 | 2 | 0 (-) | 1 |
| ai-in-package.Rmd | 8 | 6 | 2 | 5 (ran) | 0 |
| ai-quickstart.Rmd | 9 | 5 | 2 | 5 (ran) | 2 |
| B-references.Rmd | 0 | 0 | 0 | 0 (-) | 0 |
| cell-cell-communication-microenvironment.Rmd | 19 | 17 | 0 | 0 (-) | 2 |
| cell-fate-dynamics.Rmd | 26 | 9 | 14 | 0 (-) | 3 |
| cell-identity-reference-mapping.Rmd | 8 | 4 | 1 | 0 (-) | 3 |
| core-single-cell-workflow.Rmd | 18 | 18 | 0 | 0 (-) | 0 |
| data-ingress-interoperability.Rmd | 0 | 0 | 0 | 0 (-) | 0 |
| differential-analysis-interpretation.Rmd | 6 | 4 | 1 | 0 (-) | 1 |
| index.Rmd | 2 | 2 | 0 | 0 (-) | 0 |
| integration-batch-correction.Rmd | 13 | 4 | 6 | 0 (-) | 3 |
| interactive.Rmd | 3 | 1 | 0 | 0 (-) | 2 |
| milo.Rmd | 4 | 4 | 0 | 0 (-) | 0 |
| multimodal-expansion.Rmd | 16 | 2 | 13 | 0 (-) | 1 |
| part-01-foundations.Rmd | 0 | 0 | 0 | 0 (-) | 0 |
| part-02-core-workflow.Rmd | 0 | 0 | 0 | 0 (-) | 0 |
| part-03-inference.Rmd | 0 | 0 | 0 | 0 (-) | 0 |
| part-04-ai-assisted.Rmd | 0 | 0 | 0 | 0 (-) | 0 |
| phenotype-association.Rmd | 2 | 1 | 0 | 0 (-) | 1 |
| pipeline-logs.Rmd | 4 | 0 | 4 | 0 (-) | 0 |
| program-regulon-mechanism.Rmd | 15 | 11 | 0 | 0 (-) | 4 |
| spatial-context-niche.Rmd | 7 | 3 | 0 | 0 (-) | 4 |
| state-aware.Rmd | 3 | 2 | 1 | 0 (-) | 0 |
| state-priority-perturbation-sensitivity.Rmd | 13 | 6 | 7 | 1 (ran) | 0 |
| svp.Rmd | 16 | 0 | 16 | 0 (-) | 0 |
| visualization.Rmd | 24 | 22 | 0 | 0 (-) | 2 |

**Totals**: 234 chunks, 121 executed, 84 conditional, 29 off; of the 11 AI-gated chunks, 11 ran.

