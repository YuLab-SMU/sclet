# ai_report_pbmc10k_annotation

- Recorded: 2026-10-20 15:31:23.439289
- Data: 33538 genes x 11769 cells
- Ledger fingerprint at record time: `sclet-ai-3eeb268e`; now: `sclet-ai-52685d27`
- Findings: 4

## Summary

PBMC 10k 的 13 个群已完成参考注释（MonacoImmuneData / label.main）：12 个群有 >= 0.7 的主导标签，1 个群标签不一致。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| 在参考 MonacoImmuneData（label.main）下，13 个群里有 12 个群的单一标签占比 >= 0.7；cluster 9 的主导标签只占 0.22，不能给单一类型。 | `ev:annotation_monaco` | associated | ok |
| 组成（按群）：CD4+ T cells (clusters 1, 7, 9)；Monocytes (clusters 10, 12, 3)；Progenitors (clusters 11, 13)；Dendritic cells (clusters 2, 8)；NK cells (clusters 4)；B cells (clusters 5)；T cells (clusters 6)。 | `ev:annotation_monaco`, `ev:findallmarkers` | associated | ok |
| 单核细胞与 CD4+ T 细胞群各自成组（共 12 个高一致性群），说明 marker 结构与参考标签在该粒度上不矛盾。 | `ev:annotation_monaco`, `ev:findallmarkers` | associated | ok |
| cluster 9 的低一致性可能反映过渡态、双细胞或参考覆盖不足，本次无法区分。 | - | hypothesis | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:annotation_monaco` | consistent_with | `monaco` | 6c1363bc0e1f |
| `ev:findallmarkers` | associated | `markers_by_cluster` | 077cd8c55792 |

## Claim ceiling

- Ceiling: `associated` (2 independent line(s) of support)
- Reason: 2 independent lines of support; the ceiling is additionally bounded by the weakest supporting node

## Warnings and recommendations

- Warning: 参考与标签列是用户通过澄清通道记录的选择（MonacoImmuneData, label.main），不是引擎猜的。
- Warning: SingleR 平均得分 0.097 偏低；标签应视为参考映射，而非细胞身份的确认。
- Warning: 本次没有用 pruned.labels，也没有按得分阈值过滤细胞。
- Recommendation: 若要区分混合群，先跑 doublet 检测，再用细粒度参考复核。
- Recommendation: 若要给条件效应，用 run_de_test(group=, covariates=) 在配对 donor 上做。

## Limitations stated in the findings

- **在参考 MonacoImmuneData（label.main）下，13 个群里有 12 个群的单一标签占比 >= 0.7；cluster 9 的主导标签只占 0.22，不能给单一类型。**: SingleR 平均得分 0.097（0-1），分数偏低，标签置信度应按此打折。; 标签来自 label.main 粒度，比 label.fine 粗。
- **组成（按群）：CD4+ T cells (clusters 1, 7, 9)；Monocytes (clusters 10, 12, 3)；Progenitors (clusters 11, 13)；Dendritic cells (clusters 2, 8)；NK cells (clusters 4)；B cells (clusters 5)；T cells (clusters 6)。**: 组成比例是细胞数占比，不是条件间变化。
- **单核细胞与 CD4+ T 细胞群各自成组（共 12 个高一致性群），说明 marker 结构与参考标签在该粒度上不矛盾。**: 两条证据线共享同一表达矩阵：按账本的独立线规则记为 2 条，但它们不是两次独立实验——注释来自外部参考，marker 检验来自群标签。
- **cluster 9 的低一致性可能反映过渡态、双细胞或参考覆盖不足，本次无法区分。**: 没有 doublet 证据，也没有更细粒度的参考来区分过渡态与混合群。

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

### annotation `monaco`

- Method: SingleR

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

### ai_execution_plan_m5_annotation_20261020152726_8hoqtmf3_20261020152736

| step | action | status | seconds |
|---|---|---|---|
| a1 | run_annotation | completed | 9.24 |


## Verification

- Citations resolve: PASS of 2
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
