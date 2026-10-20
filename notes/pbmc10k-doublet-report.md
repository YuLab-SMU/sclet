# ai_report_pbmc10k_doublet

- Recorded: 2026-10-20 17:13:37.614152
- Data: 33538 genes x 11769 cells
- Ledger fingerprint at record time: `sclet-ai-14ceb4ad`; now: `sclet-ai-156a1131`
- Findings: 4

## Summary

PBMC 10k 的双细胞与稀有群检查：全局双细胞率 7.33%；cluster 12 是双细胞群（98.3%），cluster 9 的混合标签与双细胞无关（0%）；未发现稀有群。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| scDblFinder 在 11769 个细胞里判出 863 个 doublet（7.33%）；这是单一分类器的估计，没有第二类证据佐证，也没有任何细胞被移除。 | `ev:doublet_doublet_finder` | associated | ok |
| cluster 12 几乎全是 doublet：291/296（98.3%），远高于全局 7.33%。它的 Monaco 标签是 Monocytes，但这个标签不应被读成细胞类型——这是双细胞成团的典型形状。 | `ev:doublet_doublet_finder`, `ev:annotation_monaco`, `ev:findallmarkers` | associated | ok |
| cluster 9 的 doublet 为 0%（0/719），所以它的标签混杂（主导 CD4+ T cells 仅 22%）不是双细胞造成的。 | `ev:doublet_doublet_finder`, `ev:annotation_monaco` | associated | ok |
| 密度法在 2,350 个密度极大值上解散后得到 0 个稀有群，因此本次不主张任何稀有群。 | - | hypothesis | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:doublet_doublet_finder` | associated | - | dc0bf290ae75 |
| `ev:annotation_monaco` | consistent_with | `monaco` | 6c1363bc0e1f |
| `ev:findallmarkers` | associated | `markers_by_cluster` | 077cd8c55792 |

## Claim ceiling

- Ceiling: `associated` (3 independent line(s) of support)
- Reason: 3 independent lines of support; the ceiling is additionally bounded by the weakest supporting node

## Warnings and recommendations

- Warning: 双细胞只做了标记，没有移除任何细胞；下游统计仍包含被标为 doublet 的细胞。
- Warning: 稀有群检测得到 0 个群（q = 0.25），按设计不写证据节点，故不可引用。
- Warning: cluster 12 的 Monaco 标签（单核）来自被双细胞主导的群，报告组成时应剔除或标注。
- Recommendation: 若要得到干净组成，按 doublet 标签过滤后重跑聚类与注释。
- Recommendation: 若要主张稀有群，换更细的分辨率或更小的稀有阈值重跑，并确认有 markers/doublet 信号可分级。

## Limitations stated in the findings

- **scDblFinder 在 11769 个细胞里判出 863 个 doublet（7.33%）；这是单一分类器的估计，没有第二类证据佐证，也没有任何细胞被移除。**: 双细胞率是模型判断，不是对 doubleness 的直接测量。
- **cluster 12 几乎全是 doublet：291/296（98.3%），远高于全局 7.33%。它的 Monaco 标签是 Monocytes，但这个标签不应被读成细胞类型——这是双细胞成团的典型形状。**: 这一条修正了先前组成行里把该群计入单核细胞的做法。
- **cluster 9 的 doublet 为 0%（0/719），所以它的标签混杂（主导 CD4+ T cells 仅 22%）不是双细胞造成的。**: 排除了双细胞解释，仍不能区分过渡态与参考覆盖不足。
- **密度法在 2,350 个密度极大值上解散后得到 0 个稀有群，因此本次不主张任何稀有群。**: 零个稀有群时该动作按设计不注册证据节点，所以这句话不构成可引用的结论。

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

### rare_cells `rareq`

- Method: density
- reduction: PCA
- n_clusters: 9

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

### ai_execution_plan_m7_doublet_20261020170013_f6i7a4p3_20261020170803

| step | action | status | seconds |
|---|---|---|---|
| d1 | run_doublet_detection | completed | 466.82 |

### ai_execution_plan_m7_rare_20261020170803_ecyzkkht_20261020170808

| step | action | status | seconds |
|---|---|---|---|
| r1 | run_rare_cell_detection | completed | 4.07 |


## Verification

- Citations resolve: PASS of 3
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
