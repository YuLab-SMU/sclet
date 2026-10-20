# ai_report_kang_m2

- Recorded: 2026-10-20 14:58:17.909384
- Data: 32738 genes x 43095 cells
- Ledger fingerprint at record time: `sclet-ai-65cef876`; now: `sclet-ai-442d8ed2`
- Findings: 4

## Summary

Kang IFN-beta 的批次校正示范闭环已完成：设计确认 → fastMNN/Harmony 两条路线 → 路线比较（发现 trade-off）→ 配对条件对比。全部结论的最高证据等级是 associated。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| 两个批次在未校正的 PCA 上有明显分离（batch_mixing 0.9081）；fastMNN 与 Harmony 都把它推高到 0.9734 / 0.9913，均在校正自身的嵌入上测得。 | `ev:integration_comparison` | associated | ok |
| 校正存在 trade-off：batch_mixing 最优的路线（route_3）把 biological_preservation 从基线的 0.0751 降到 0.0598，而 fastMNN 在两项上同时改善（0.0842），代价是 59.5 s 对 9.4 s。 | `ev:integration_comparison` | associated | ok |
| stim vs ctrl 的配对差异检验（Pseudobulk limma (condition contrast, paired by donor)，16 个 pseudobulk 样本，设计 ~ condition + donor）得到 6537 行结果，其中 4254 行 p_adj < 0.05、1074 行在 stim 上调。 | `ev:detest_condition_stim_vs_ctrl` | associated | ok |
| stim 下上调的基因可能反映干扰素刺激反应；需要通路富集（run_enrichment）或基于参考的注释才能支持这一读法。 | - | hypothesis | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:integration_comparison` | associated | `raw`, `fastmnn`, `harmony` | e436fae9a091 |
| `ev:detest_condition_stim_vs_ctrl` | associated | `stim_vs_ctrl_paired` | e12381600a4d |

## Claim ceiling

- Ceiling: `associated` (2 independent line(s) of support)
- Reason: 2 independent lines of support; the ceiling is additionally bounded by the weakest supporting node

## Warnings and recommendations

- Warning: stim 只在 batch2 中存在，batch 与 condition 部分混杂：任何'条件效应'都受批次残留影响。
- Warning: biological_preservation 取自 condition 列的 eta^2，不是细胞类型保真度。
- Warning: cluster_stability 指标在本版本未实现（stability_resampling_skipped_in_phase_c）。
- Recommendation: 若要判断哪条路线更适合下游分析，可在两条路线上各跑一遍 find_clusters + run_de_test，比较 marker 的一致性。
- Recommendation: 若要给 IFN-beta 响应一个生物学读法，跑 run_enrichment（db 由用户指定）。

## Limitations stated in the findings

- **两个批次在未校正的 PCA 上有明显分离（batch_mixing 0.9081）；fastMNN 与 Harmony 都把它推高到 0.9734 / 0.9913，均在校正自身的嵌入上测得。**: batch_mixing = 1 - eta^2，eta^2 取前 10 维的均值；它是结构描述，不是校正优劣的结论。
- **校正存在 trade-off：batch_mixing 最优的路线（route_3）把 biological_preservation 从基线的 0.0751 降到 0.0598，而 fastMNN 在两项上同时改善（0.0842），代价是 59.5 s 对 9.4 s。**: condition 与 batch 部分混杂（stim 只存在于 batch2），所以两个指标不独立。; cluster_stability 指标尚未实现（status = not_available）。
- **stim vs ctrl 的配对差异检验（Pseudobulk limma (condition contrast, paired by donor)，16 个 pseudobulk 样本，设计 ~ condition + donor）得到 6537 行结果，其中 4254 行 p_adj < 0.05、1074 行在 stim 上调。**: 这是配对差异的统计检验，不说明机制；行数按记录的阈值（min.pct 0.01，logFC >= 0.1）统计。
- **stim 下上调的基因可能反映干扰素刺激反应；需要通路富集（run_enrichment）或基于参考的注释才能支持这一读法。**: 本次未做富集或注释，无法把上调基因映射到具体生物学过程。

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


## Verification

- Citations resolve: PASS of 2
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
