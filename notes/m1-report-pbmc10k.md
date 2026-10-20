# ai_report_pbmc10k_m1

- Recorded: 2026-10-20 14:57:10.540367
- Data: 33538 genes x 11769 cells
- Ledger fingerprint at record time: `sclet-ai-4d3dfbea`; now: `sclet-ai-59268ce2`
- Findings: 4

## Summary

PBMC 10k 的两轮分析（建图/聚类/标记基因 + GO/KEGG 注释）已完成，所有结论的最高证据等级是 associated：marker 与由它派生的富集属于同一条证据线。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| 在 resolution 0.8、seed 1 下得到 13 个 louvain 群；对全部群做标记基因检验，共 12588 条基因-群记录（min.pct 0.1，logFC >= 0.25，最小群 39 个细胞）。 | `ev:findallmarkers` | associated | ok |
| marker 基因在 GO BP 上得到 6430 条富集条目，覆盖 13 个群。 | `ev:enrichment_enrich_go` | associated | ok |
| 同一批 marker 转成 ENTREZID 后在 KEGG 上得到 645 条富集条目，覆盖 13 个群。 | `ev:enrichment_enrich_kegg` | associated | ok |
| 这些群可能对应不同的 PBMC 细胞类型（例如含一个约 39 个细胞的小群），但本次没有做基于参考的注释。 | - | hypothesis | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:findallmarkers` | associated | `markers_by_cluster` | 077cd8c55792 |
| `ev:enrichment_enrich_go` | consistent_with | `ev:findallmarkers` | b355e1f4c98e |
| `ev:enrichment_enrich_kegg` | consistent_with | `ev:findallmarkers` | b355e1f4c98e |

## Claim ceiling

- Ceiling: `associated` (1 independent line(s) of support)
- Reason: all supporting evidence resolves to 1 independent line(s) of support; a claim cannot exceed 'associated'

## Warnings and recommendations

- Warning: 只有一个样本，没有条件对比，也没有 donor/batch 结构，因此不能做条件或批次层结论。
- Warning: 没有做基于参考的细胞类型注释（需要用户指定 celldex 参考）。
- Recommendation: 若需要细胞类型，用用户指定的 celldex 参考跑 run_annotation。
- Recommendation: 若要评价条件效应，先确认 design 角色，再用 run_de_test(group=, covariates=)。

## Limitations stated in the findings

- **在 resolution 0.8、seed 1 下得到 13 个 louvain 群；对全部群做标记基因检验，共 12588 条基因-群记录（min.pct 0.1，logFC >= 0.25，最小群 39 个细胞）。**: 只说明各群相对其余细胞的表达差异，不说明群里是什么细胞。
- **marker 基因在 GO BP 上得到 6430 条富集条目，覆盖 13 个群。**: 条目名与基因列表留在 enrichment 状态记录里，证据节点按设计不含名称。
- **同一批 marker 转成 ENTREZID 后在 KEGG 上得到 645 条富集条目，覆盖 13 个群。**: KEGG 与 GO 是同一份 marker 表的两种注释，不是两次独立测量。
- **这些群可能对应不同的 PBMC 细胞类型（例如含一个约 39 个细胞的小群），但本次没有做基于参考的注释。**: 把 cluster_N 读成细胞类型是本项目的禁区之一。; 没有注释证据；群身份只能由 marker 推测。

## Recorded analyses

### clustering `louvain_clusters`

- Method: igraph::cluster_louvain
- reduction: PCA
- n_clusters: 13
- resolution: 0.8
- seed: 1

### detest `markers_by_cluster`

- Method: FindAllMarkers
- n_markers: 12588
- min.pct: 0.1
- logfc.threshold: 0.25

### detest `markers_entrez`

- Method: clusterProfiler::bitr (derived from markers_by_cluster)
- n_markers: 12213

### enrichment `enrich_go`

- Method: clusterProfiler::compareCluster
- db: GO
- n_enriched: 6430

### enrichment `enrich_kegg`

- Method: clusterProfiler::compareCluster
- db: KEGG
- n_enriched: 645

## Executed steps (from the execution records)

### ai_execution_plan_m1_pipeline_20261020133318_83dkind4_20261020133400

| step | action | status | seconds |
|---|---|---|---|
| s1 | check_qc | completed | 0.09 |
| s2 | normalize_data | completed | 2.16 |
| s3 | find_variable_features | completed | 0.33 |
| s4 | scale_data | completed | 8.9 |
| s5 | run_pca | completed | 26.37 |
| s6 | find_neighbors | completed | 1.76 |
| s7 | find_clusters | completed | 1.36 |

### ai_execution_plan_m1_markers_20261020133400_f4b51qir_20261020133403

| step | action | status | seconds |
|---|---|---|---|
| m1 | run_de_test | completed | 3.1 |

### ai_execution_plan_m1_enrichment_go_20261020133405_dreqobj9_20261020133513

| step | action | status | seconds |
|---|---|---|---|
| g1 | run_enrichment | completed | 67.35 |

### ai_execution_plan_m1_enrichment_kegg_20261020133513_vehtih7h_20261020133634

| step | action | status | seconds |
|---|---|---|---|
| k1 | run_enrichment | completed | 80.24 |


## Verification

- Citations resolve: PASS of 3
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
