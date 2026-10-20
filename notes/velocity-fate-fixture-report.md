# ai_report_velocity_fate_fixture

- Recorded: 2026-10-20 19:05:22.628902
- Data: 200 genes x 120 cells
- Ledger fingerprint at record time: `sclet-ai-0cf09045`; now: `sclet-ai-261f8a7a`
- Findings: 2

## Summary

合成 velocity fixture（120 细胞，spliced/unspliced/counts）上跑通 scvelo + CellRank；fate 退化为单一终态（1 个 terminal state，全部细胞概率 1），所以它证明链路可用，不产生生物学结论。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| scvelo 估计出速度（mode_code 1），写入 4 个 colData 字段（pseudotime/confidence/root/end）。 | `ev:velocity_velo` | associated | ok |
| CellRank 只找到 1 个 terminal state，全部细胞吸收概率为 1：fate 结果没有区分度。 | `ev:fate_fate` | associated | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:velocity_velo` | consistent_with | - | dc0bf290ae75 |
| `ev:fate_fate` | consistent_with | - | dc0bf290ae75 |

## Claim ceiling

- Ceiling: `associated` (1 independent line(s) of support)
- Reason: all supporting evidence resolves to 1 independent line(s) of support; a claim cannot exceed 'associated'

## Warnings and recommendations

- Warning: 整段链条跑在合成 fixture 上：真实 velocity 导出（hermann sperm、regvelo smoke）都没有 counts assay。
- Warning: fate 退化说明“链路通”不等于“结果有意义”；要得到可解释的 fate 需要多终态的轨迹与更细的群。
- Recommendation: 若要真实 fate，先在带 counts 的真实数据上跑 velocity，再用多终态的轨迹。

## Limitations stated in the findings

- **scvelo 估计出速度（mode_code 1），写入 4 个 colData 字段（pseudotime/confidence/root/end）。**: 速度是模型估计；本次没有写出 velocity reduction（n_reductions_written = 0）。
- **CellRank 只找到 1 个 terminal state，全部细胞吸收概率为 1：fate 结果没有区分度。**: 数据是合成的，且 stage 三分后只有一个 macrostate；这条不能读成生物学结论。; CellRank 警告只有一个 lineage，改用 stationary distribution。

## Recorded analyses

### trajectory `fate`

- Method: CellRank
- reduction: PCA

### velocity `velo`

- Method: velociraptor::scvelo
- reduction: UMAP

## Executed steps (from the execution records)

### ai_execution_plan_m9_f_prep_20261020190221_dh0zctb2_20261020190225

| step | action | status | seconds |
|---|---|---|---|
| p1 | run_pca | completed | 0.14 |
| p2 | run_umap | completed | 2.91 |

### ai_execution_plan_m9_f_velocity_20261020190225_6on1o0nr_20261020190234

| step | action | status | seconds |
|---|---|---|---|
| v1 | run_velocity | completed | 8 |

### ai_execution_plan_m9_f_fate_20261020190234_90x4jotu_20261020190317

| step | action | status | seconds |
|---|---|---|---|
| f1 | run_fate_analysis | completed | 42.75 |


## Verification

- Citations resolve: PASS of 2
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
