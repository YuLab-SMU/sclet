# ai_report_kang_milo

- Recorded: 2026-10-20 16:20:10.018497
- Data: 32738 genes x 43095 cells
- Ledger fingerprint at record time: `sclet-ai-5fd55a43`; now: `sclet-ai-1585ff59`
- Findings: 4

## Summary

Kang IFN-beta 的差异丰度分析完成：3414 个邻域中 1821 个在 stim vs ctrl 下显著（FDR<0.05）；置换标签的对照为 0/3414。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| stim 与 ctrl 的邻域丰度比较（~0 + condition + batch，对比 stim - ctrl）在 3414 个邻域里有 1821 个 FDR < 0.05（中位 p 0.0152）。 | `ev:milo_da_stim` | associated | ok |
| 同一批细胞的配对 pseudobulk 差异表达（~ condition + donor，16 个样本）有 4254 行 p_adj < 0.05；组成与表达两个层面都检出 stim 效应。 | `ev:detest_condition_stim_vs_ctrl`, `ev:milo_da_stim` | associated | ok |
| 把 4 个配对 donor 的 condition 标签对调后重跑，0/3414 个邻域显著（中位 p 0.7586），说明上面的信号不是标签顺序造成的。 | `ev:milo_da_stim` | associated | ok |
| 邻域丰度的变化可能反映单核/树突等群体的相对扩张或收缩，本次无法区分扩张与状态转换。 | - | hypothesis | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:milo_da_stim` | associated | - | milo_a55e260 |
| `ev:detest_condition_stim_vs_ctrl` | associated | `stim_vs_ctrl_paired` | e12381600a4d |

## Claim ceiling

- Ceiling: `associated` (2 independent line(s) of support)
- Reason: 2 independent lines of support; the ceiling is additionally bounded by the weakest supporting node

## Warnings and recommendations

- Warning: 配对 donor 结构未进入 Milo 设计（formula ~0 + condition + batch），因此这是未配对的邻域检验；配对版本需要 GLMM + donor 随机效应。
- Warning: 样本里有一个单细胞样本（batch1:unknown），已随设计进入 TMM 归一化。
- Warning: 差异丰度描述组成变化，不描述表达变化；两者不能互相替代。
- Recommendation: 若要利用配对结构，用 glmm = TRUE + random_effect = 'donor' 重跑并比较显著邻域数。
- Recommendation: 若要解释邻域变化，把 nhood groups 映射到注释标签后再看方向。

## Limitations stated in the findings

- **stim 与 ctrl 的邻域丰度比较（~0 + condition + batch，对比 stim - ctrl）在 3414 个邻域里有 1821 个 FDR < 0.05（中位 p 0.0152）。**: 邻域丰度变化不等于表达变化；它描述的是细胞组成的位移。; 配对 donor 结构没有建模（Milo 设计为 ~ 0 + condition + batch，未用 donor 随机效应）。
- **同一批细胞的配对 pseudobulk 差异表达（~ condition + donor，16 个样本）有 4254 行 p_adj < 0.05；组成与表达两个层面都检出 stim 效应。**: 两个证据来自同一表达矩阵；按账本规则记为 2 条独立线，但不是两次独立实验。
- **把 4 个配对 donor 的 condition 标签对调后重跑，0/3414 个邻域显著（中位 p 0.7586），说明上面的信号不是标签顺序造成的。**: 置换是在 donor 内做的，只打乱条件标签，保留批次与细胞组成。
- **邻域丰度的变化可能反映单核/树突等群体的相对扩张或收缩，本次无法区分扩张与状态转换。**: 没有做邻域层面的细胞类型分解（Milo nhood groups 未纳入），也没有 doublet 证据。

## Recorded analyses

### clustering `louvain_clusters`

- Method: igraph::cluster_louvain
- reduction: PCA
- n_clusters: 9
- resolution: 0.5
- seed: 1

### detest `stim_vs_ctrl_paired`

- Method: Pseudobulk limma (condition contrast, paired by donor)
- design: ~ condition + donor
- adjustment: pseudobulk per donor and limma, design ~ condition + donor
- n_pseudobulk_samples: 16
- n_markers: 6537

### integration `fastmnn`

- Method: batchelor::batchCorrect

### integration `harmony`

- Method: harmony::RunHarmony
- reduction: corrected
- batch: batch

### integration `integration_comparison`

- Method: sclet_ai_route_comparison

| route | batch mixing (higher better) | biological preservation (higher better) | seconds | embedding | baseline |
|---|---|---|---|---|---|
| raw | 0.9081 | 0.0751 | 0 | PCA | yes |
| fastmnn | 0.9734 | 0.0842 | 59.45 | corrected | no |
| harmony | 0.9913 | 0.0598 | 9.42 | HARMONY | no |

- Batch mixing improved: TRUE; best route kept the biological signal: FALSE; trade-off: TRUE

### milo `da_stim`

- Method: miloR
- adjustment: BH
- reduction: PCA
- seed: 2025

## Executed steps (from the execution records)

### ai_execution_plan_m2_design_20261020091754_94a4dcup_20261020091755

| step | action | status | seconds |
|---|---|---|---|
| d1 | confirm_design | completed | 0.08 |

### ai_execution_plan_m2_route_fastmnn_20261020091755_idtpofgl_20261020091855

| step | action | status | seconds |
|---|---|---|---|
| i1 | run_integration | completed | 59.6 |

### ai_execution_plan_m2_route_harmony_20261020091855_qpi6halj_20261020091905

| step | action | status | seconds |
|---|---|---|---|
| i2 | run_integration | completed | 9.56 |

### ai_execution_plan_m2_compare_20261020091905_fqycjy1q_20261020091906

| step | action | status | seconds |
|---|---|---|---|
| c1 | compare_integration_routes | completed | 0.28 |

### ai_execution_plan_m2_condition_de_20261020091907_w0jo70q7_20261020091909

| step | action | status | seconds |
|---|---|---|---|
| p1 | run_de_test | completed | 2.03 |

### ai_execution_plan_m6_design_20261020160632_x0wjjkh5_20261020160633

| step | action | status | seconds |
|---|---|---|---|
| d1 | confirm_design | completed | 0.08 |

### ai_execution_plan_m6_milo_20261020160633_upetk5w4_20261020160823

| step | action | status | seconds |
|---|---|---|---|
| m1 | run_milo | completed | 108.49 |


## Verification

- Citations resolve: PASS of 2
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
