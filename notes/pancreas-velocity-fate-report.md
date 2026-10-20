# ai_report_pancreas_velocity_fate

- Recorded: 2026-10-20 23:26:53.129446
- Data: 27998 genes x 3696 cells
- Ledger fingerprint at record time: `sclet-ai-4d643158`; now: `sclet-ai-54e8042d`
- Findings: 4

## Summary

真实数据上的速度-轨迹-命运链跑通（3696 细胞）：5 条 lineage、3 个终态，去向分布 0.107/0.599/0.294。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| 在 scvelo pancreas 数据（3696 个细胞，真实 spliced/unspliced）上跑通 scVelo deterministic 速度估计，projected 到 UMAP，45 s。 | `ev:velocity_velo` | associated | ok |
| 从 Ngn3 low EP 出发的 slingshot 拟合出 5 条 lineage，覆盖 8 个群；pseudotime 只表达相对顺序。 | `ev:trajectory_pancreas_traj` | associated | ok |
| CellRank 估计出 3 个 terminal state；按每个细胞吸收概率取 argmax，去向分布是 0.107 / 0.599 / 0.294，三个 fate 维度的平均吸收概率是 0.214 / 0.524 / 0.262。 | `ev:fate_fate`, `ev:velocity_velo`, `ev:trajectory_pancreas_traj` | associated | ok |
| 证据里的 terminal_state_fractions（0.333/0.333/0.333）是 90 个终态细胞之间的占比，不是全体 3696 个细胞的去向；读去向要用 destination_fractions。 | `ev:fate_fate` | observed | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:velocity_velo` | consistent_with | - | dc0bf290ae75 |
| `ev:trajectory_pancreas_traj` | consistent_with | - | dc0bf290ae75 |
| `ev:fate_fate` | consistent_with | - | dc0bf290ae75 |

## Claim ceiling

- Ceiling: `associated` (1 independent line(s) of support)
- Reason: all supporting evidence resolves to 1 independent line(s) of support; a claim cannot exceed 'associated'

## Warnings and recommendations

- Warning: terminal_state_fractions 只是终态细胞之间的占比（90/3696 个细胞），按 1/k 构造；去向分布另存 destination_fractions。
- Warning: 终态与群的编码是匿名的，报告没有给出 alpha/beta/delta 的名字。
- Recommendation: 若要把终态映射到细胞类型，用记录里的 group_codes 顺序对照注释列。
- Recommendation: 若要判断命运是否可信，用另一套速度参数（stochastic/dynamical）复跑并比较去向分布。

## Limitations stated in the findings

- **在 scvelo pancreas 数据（3696 个细胞，真实 spliced/unspliced）上跑通 scVelo deterministic 速度估计，projected 到 UMAP，45 s。**: 速度是模型估计；counts 取自 spliced 计数（velocity 导出的惯例），不是独立测序。
- **从 Ngn3 low EP 出发的 slingshot 拟合出 5 条 lineage，覆盖 8 个群；pseudotime 只表达相对顺序。**: 起点取自数据集自带注释（Ngn3 low EP 是内分泌谱系的已知起点），不是拟合推出。
- **CellRank 估计出 3 个 terminal state；按每个细胞吸收概率取 argmax，去向分布是 0.107 / 0.599 / 0.294，三个 fate 维度的平均吸收概率是 0.214 / 0.524 / 0.262。**: 终态用 cluster_1..3 的匿名编码表示，本次没有把编码映射回 alpha/beta/delta。; 吸收概率是模型估计，不是命运观测。
- **证据里的 terminal_state_fractions（0.333/0.333/0.333）是 90 个终态细胞之间的占比，不是全体 3696 个细胞的去向；读去向要用 destination_fractions。**: 这是字段语义的读数，不是生物学结论。

## Recorded analyses

### trajectory `pancreas_traj`

- Method: slingshot
- reduction: UMAP
- seed: 2025

### trajectory `fate`

- Method: CellRank
- reduction: PCA

### velocity `velo`

- Method: velociraptor::scvelo
- reduction: UMAP

## Executed steps (from the execution records)

### ai_execution_plan_m11_prep_20261020224527_0v58uaff_20261020224540

| step | action | status | seconds |
|---|---|---|---|
| p1 | run_pca | completed | 5.77 |
| p2 | run_umap | completed | 6.24 |

### ai_execution_plan_m11_velocity_20261020224540_gia79v2k_20261020224625

| step | action | status | seconds |
|---|---|---|---|
| v1 | run_velocity | completed | 44.78 |

### ai_execution_plan_m11_trajectory_20261020225847_qe1w65xz_20261020225854

| step | action | status | seconds |
|---|---|---|---|
| t1 | run_trajectory | completed | 5.73 |

### ai_execution_plan_m11_fate_20261020230257_05q4ypbi_20261020231157

| step | action | status | seconds |
|---|---|---|---|
| f1 | run_fate_analysis | completed | 538.97 |


## Verification

- Citations resolve: PASS of 3
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
