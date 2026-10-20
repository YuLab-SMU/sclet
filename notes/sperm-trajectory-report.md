# ai_report_sperm_trajectory

- Recorded: 2026-10-20 19:05:21.399277
- Data: 2000 genes x 300 cells
- Ledger fingerprint at record time: `sclet-ai-10395ceb`; now: `sclet-ai-4ba11068`
- Findings: 3

## Summary

Hermann sperm 数据（300 细胞，5 群）上 slingshot 得到 1 条 lineage，起点 cluster_1（Pachytene spermatocytes，数据里最早的分化阶段）。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| 从 cluster_1 出发，slingshot 拟合出 1 条 lineage，覆盖 5 个群；pseudotime 只用相对顺序，不是绝对时间。 | `ev:trajectory_sperm_traj` | associated | ok |
| 第一名 lineage 的 pseudotime 中位数 18.1796、IQR 9.2671，说明细胞沿该 lineage 分布较宽，不是一条细线。 | `ev:trajectory_sperm_traj` | associated | ok |
| 这些群的先后顺序可能对应精子发生过程（Pachytene → round spermatids），但本次没有用速度或外部标记验证方向。 | - | hypothesis | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:trajectory_sperm_traj` | consistent_with | - | dc0bf290ae75 |

## Claim ceiling

- Ceiling: `associated` (1 independent line(s) of support)
- Reason: all supporting evidence resolves to 1 independent line(s) of support; a claim cannot exceed 'associated'

## Warnings and recommendations

- Warning: RNA velocity 在本数据上没有跑成：velociraptor 需要 counts assay，而这份对象只有 spliced/unspliced/logcounts；现在计划层会直接拒绝（counts_assay_missing）而不是等到 Python 里崩。
- Warning: 因此 slingshot 的方向性只是群顺序，没有速度证据支持。
- Recommendation: 若要方向证据，提供带 counts 的 spliced/unspliced 对象后重跑 run_velocity。

## Limitations stated in the findings

- **从 cluster_1 出发，slingshot 拟合出 1 条 lineage，覆盖 5 个群；pseudotime 只用相对顺序，不是绝对时间。**: 起点是建模选择：cluster_1 是数据里最早阶段，不是推断出的根。; 只有 5 个群、300 个细胞，曲线形状对群划分敏感。
- **第一名 lineage 的 pseudotime 中位数 18.1796、IQR 9.2671，说明细胞沿该 lineage 分布较宽，不是一条细线。**: 这些是归一化后的相对量，不能跨数据集比较。
- **这些群的先后顺序可能对应精子发生过程（Pachytene → round spermatids），但本次没有用速度或外部标记验证方向。**: velocity 因缺少 counts assay 未能运行（见警告），方向只有群顺序这一条线索。

## Recorded analyses

### clustering `louvain_clusters`

- Method: igraph::cluster_louvain
- reduction: PCA
- n_clusters: 5
- resolution: 0.5

### trajectory `sperm_traj`

- Method: slingshot
- reduction: UMAP
- seed: 2025

## Executed steps (from the execution records)

### ai_execution_plan_m9_inspect_ledger_20261020185956_l4n8gij4_20261020185956

| step | action | status | seconds |
|---|---|---|---|
| l1 | inspect_ledger | completed | 0.06 |

### ai_execution_plan_m9_prep_20261020185956_2szmprhq_20261020190001

| step | action | status | seconds |
|---|---|---|---|
| p1 | run_pca | completed | 0.96 |
| p2 | run_umap | completed | 2.88 |

### ai_execution_plan_m9_velocity_20261020190001_4dysmp9k_20261020190001

| step | action | status | seconds |
|---|---|---|---|
| v1 | run_velocity | failed | 0.02 |

### ai_execution_plan_m9_trajectory_20261020190157_s2rdklmy_20261020190159

| step | action | status | seconds |
|---|---|---|---|
| t1 | run_trajectory | completed | 1.32 |


## Verification

- Citations resolve: PASS of 1
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
