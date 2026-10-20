# Kang IFN-beta 批次校正报告（M2 验收物：示范闭环）

> 数据：Kang et al. 2018 IFN-beta PBMC（GSE96583），32738 基因 × 43095 细胞；两个测序批次、8 个已配对 donor（stim 仅存在于 batch2）。
> 批次规模：batch1 14030 细胞（全 ctrl），batch2 29065 细胞（ctrl 14932 + stim 14133）。
> 路径：基线 pipeline → 设计确认 → fastMNN → Harmony → 路线比较 → 配对条件对比，五轮分轮计划，每轮重新校验。
> 账本指纹（记录本报告之前）：`sclet-ai-7b37cb3c`；记录后 `sclet-ai-17c7bf0e`。
> 报告记录 id：`ai_report_kang_m2`；sclet 1.0.2 @ `9400e1d`（+ M2 路线度量改动）；生成时间 2026-10-20 10:37。

## 路线比较（每个数字都在路线自身的嵌入上测得）

| 路线 | batch_mixing ↑ | biological_preservation ↑ | runtime (s) ↓ | 嵌入 |
|---|---|---|---|---|
| raw | 0.9081 | 0.0751 | 0 | PCA |
| fastmnn | 0.9734 | 0.0842 | 59.4 | corrected |
| harmony | 0.9913 | 0.0598 | 9.4 | HARMONY |

比较的结论（记录于 `integration_comparison` 记录与 `ev:integration_comparison` 节点）：batch mixing 改善 = **TRUE**，最优路线的 biological preservation 保留 = **FALSE**，trade-off = **TRUE**。

注：batch_mixing = 1 − η²（前 10 维，对 batch 列），biological_preservation = η²（对 condition 列）；
η² 衡量嵌入中该列能解释的方差比例。基线路线（raw）在对象自己的 PCA 上测得。

## 配对条件对比（stim vs ctrl）

- 方法：Pseudobulk limma (condition contrast, paired by donor)（`run_de_test(group=, covariates=)` 第一次在真实数据上执行）。
- 设计：~ condition + donor；pseudobulk 样本 16 个（8 donor × 2 条件）。
- 结果行数 6537（按记录阈值 min.pct 0.01、logFC ≥ 0.1 统计）。

## 结论（等级由引擎按证据判定）

| 结论 | 证据引用 | 记录等级 | 引擎审计 |
|---|---|---|---|
| 两个批次在未校正的 PCA 上有明显分离（batch_mixing 0.9081）；fastMNN 与 Harmony 都把它推高到 0.9734 / 0.9913，均在校正自身的嵌入上测得。 | `ev:integration_comparison` | associated | engine: ok |
| 校正存在 trade-off：batch_mixing 最优的路线（route_3）把 biological_preservation 从基线的 0.0751 降到 0.0598，而 fastMNN 在两项上同时改善（0.0842），代价是 59.5 s 对 9.4 s。 | `ev:integration_comparison` | associated | engine: ok |
| stim vs ctrl 的配对差异检验（Pseudobulk limma (condition contrast, paired by donor)，16 个 pseudobulk 样本，设计 ~ condition + donor）得到 6537 行结果，其中 4254 行 p_adj < 0.05、1074 行在 stim 上调。 | `ev:detest_condition_stim_vs_ctrl` | associated | engine: ok |
| stim 下上调的基因可能反映干扰素刺激反应；需要通路富集（run_enrichment）或基于参考的注释才能支持这一读法。 | — | hypothesis | engine: ok |

### 为什么最高只到 `associated`

`check_claim_ceiling()`：`associated`（2 independent lines of support; the ceiling is additionally bounded by the weakest supporting node）。

两条证据线（路线比较、配对 DE）互相独立，但每个节点自身只有一条支撑，且都被标注为
`associated`——ceiling 取最弱的一环。要到 `consistent_with`，需要一条能承载它的独立证据。

## scVI：已跑通，但不在上表（细胞集不同）

roadmap 的第三条路线 scVI 本次也执行了：basilisk env 首次构建 + 在 **9000 细胞子集**（每 batch×condition 抽样）上训练，batch_mixing **0.989**、biological_preservation **0.06111**、墙钟 2877.99 s（含 env 构建）。

它没有并入上表：上表三条路线都在全部 43095 个细胞上测得，而 scVI 的数字来自另一个细胞集（9000 个），
按 sclet-provenance 的规则（不同 fixture 的分析不可合并统计）只作范围说明。env 已构建完成，
全量对象上的 scVI 现在可行（CPU 上预计数小时量级）。

## 明确没有说的话

- 没有推荐某条路线「更好」：比较只报告 trade-off（一项改善、一项代价），最终取舍是用户的建模决定。
- 没有把 biological_preservation 当作细胞类型保真度：它只看 condition 列的 η²。
- 没有给 IFN-beta 响应一个生物学读法：未做富集/注释，相关表述为 `hypothesis`。
- 没有说校正消除了批次效应：batch_mixing < 1，且 condition 与 batch 部分混杂。

## 复验方式

1. `readRDS(m2-report.rds)`，只读账本重算每个数字与等级，不读本文。
2. 断言：引用的证据 id 可解析；记录等级 = 引擎重算等级；路线数字 = 比较记录里的数字；配对设计字符串可从 DE 记录重算。
复验记录见 `notes/m2-verification-kang.md`。

## 下一步（建议，未执行）

- scVI 已在 9000 细胞子集上跑通（env 已构建）；全量对象上的第三条路线可随时补跑。
- 若要判断路线对下游的影响：在各路线上跑聚类并比较 marker 一致性。

> 警告：stim 只在 batch2 中存在，batch 与 condition 部分混杂：任何'条件效应'都受批次残留影响。
> 警告：biological_preservation 取自 condition 列的 eta^2，不是细胞类型保真度。
> 警告：cluster_stability 指标在本版本未实现（stability_resampling_skipped_in_phase_c）。
