# ai_report_pbmc10k_cci

- Recorded: 2026-10-20 17:57:32.262045
- Data: 33538 genes x 11769 cells
- Ledger fingerprint at record time: `sclet-ai-5cbd8560`; now: `sclet-ai-50126a6b`
- Findings: 3

## Summary

PBMC 10k 的细胞通讯分析完成（CellChat，用户经澄清通道指定）：12 个发送群、13 个接收群、733 对显著互作。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| CellChat 在 12 个发送群与 13 个接收群之间推断出 733 对显著配体-受体互作（概率中位数 0.0134，最大 0.1207）。 | `ev:communication_cellchat` | associated | ok |
| 该结果把 cluster 12 当作发送/接收群之一，而它的双细胞比例是 98.3%，因此涉及该群的互作不应被读成生物学信号。 | `ev:communication_cellchat`, `ev:doublet_doublet_finder` | associated | ok |
| 互作可能集中在髓系与淋巴系群之间，但本次没有引用任何配体、受体或通路名称。 | - | hypothesis | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:communication_cellchat` | consistent_with | - | cci_d63653d5 |
| `ev:doublet_doublet_finder` | associated | - | dc0bf290ae75 |

## Claim ceiling

- Ceiling: `associated` (2 independent line(s) of support)
- Reason: 2 independent lines of support; the ceiling is additionally bounded by the weakest supporting node

## Warnings and recommendations

- Warning: 引擎由用户在澄清通道里记录为 CellChat；本机只装了 CellChat，所以这个选择也受环境限制。
- Warning: 证据节点不含配体/受体/通路名称（设计如此），因此不能从证据节点直接命名通路。
- Warning: cluster 12（98.3% 双细胞）参与了互作推断，结果含污染。
- Recommendation: 若要干净的通讯结果，先按 doublet 标签过滤、重聚类，再重跑 run_cci。
- Recommendation: 若要比较引擎，用 CellPhoneDB/NicheNet 在同一细胞集上重跑——它们不是同一种测量，不能互相印证。

## Limitations stated in the findings

- **CellChat 在 12 个发送群与 13 个接收群之间推断出 733 对显著配体-受体互作（概率中位数 0.0134，最大 0.1207）。**: 互作强度是 CellChat 模型的估计；换引擎（CellPhoneDB/NicheNet）会给出不同的对与数。; 733 是 CellChat 自己检验后返回的行数，不是从全部配体-受体对里筛出的比例。
- **该结果把 cluster 12 当作发送/接收群之一，而它的双细胞比例是 98.3%，因此涉及该群的互作不应被读成生物学信号。**: 这是两个记录之间的交叉核对；要得到干净结果需先过滤双细胞再重跑。
- **互作可能集中在髓系与淋巴系群之间，但本次没有引用任何配体、受体或通路名称。**: 证据节点按设计只含计数与分位数；名称在 get_cci() 的表里，读它们需要外部知识而不是记录本身。

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

### communication `cellchat`

- Method: CellChat

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

### ai_execution_plan_m8_cci_20261020175204_tlgp6mas_20261020175359

| step | action | status | seconds |
|---|---|---|---|
| c1 | run_cci | completed | 112.96 |


## Verification

- Citations resolve: PASS of 2
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
