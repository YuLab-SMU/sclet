# M12-a 可靠性记录（自动生成，勿手改）

> 生成时间：2026-10-21 13:48:09
> 由 `R/actions/classify-reliability.R` 从 JSONL 记录生成，失败形状按封闭枚举分类（见该文件）。

## 判读规则

- 分母只计 **measured** 的 trial；`not_measured`（崩溃且没写记录）两边都不算，
  否则会把「没测到」记成「失败」或「成功」。
- 分类取 **最靠后的那一层**：空计划又被拒 → `plan_empty`；跑完但报告自检不过 →
  `report_unverifiable`，不是 `none`。
- 成功率带 Wilson 95% 区间。**不报裸百分比**：trial 数少时区间宽到会让人误读点估计。

## 汇总

| 记录文件 | trial | measured | 完成 | 成功率 | 95% 区间 | commit |
|---|---|---|---|---|---|---|
| m12a-batch-5947e4b.jsonl | 12 | 12 | 11 | 91.7% | [64.6%, 98.5%] | 5947e4b |
| realdata_kang_ifnb_bench.jsonl | 42 | 42 | 28 | 66.7% | [51.6%, 79.0%] | 4400568, fc3b7ff, a213a38, 60afa7c, 98aa395, 8d88a3b, b91ac54, ef8e4e6, 2c52dbb, 04c20aa, 80caa2e, d174367, 3e4ed58, c84a899, a8b1cc7, 2bf7b4f, 15fe726, 568a1da, 5947e4b |

## 失败形状计数

| 记录文件 | none | plan_empty | plan_shape | plan_knowledge | plan_refused | execute_failed | report_unverifiable | record_error | not_measured |
|---|---|---|---|---|---|---|---|---|---|
| m12a-batch-5947e4b.jsonl | 11 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| realdata_kang_ifnb_bench.jsonl | 28 | 2 | 3 | 1 | 6 | 2 | 0 | 0 | 0 |

## 每一条失败的具体形状

### m12a-batch-5947e4b.jsonl（1 条失败）

| # | failure_class | 层 | 记录到的错误 |
|---|---|---|---|
| 1 | plan_knowledge | plan | prerequisites not met for action round2_step_1: contrast_invalid: ident.1/ident.2 must be levels of 'condition': ctrl, stim |

### realdata_kang_ifnb_bench.jsonl（14 条失败）

| # | failure_class | 层 | 记录到的错误 |
|---|---|---|---|
| 1 | plan_refused | plan | invalid params for action not_planned_1: missing required parameter(s): ref, labels \| invalid params for action not_planned_2: missing required parameter(s): ba |
| 2 | plan_empty | plan | plan_empty: the plan has no actions. A plan with zero steps performs no analysis and cannot be executed; re-propose the plan, or if the goal genuinely needs no  |
| 3 | execute_failed | execute | execution reported failed |
| 4 | execute_failed | execute | execution reported failed |
| 5 | plan_refused | plan | prerequisites not met for action step_6: design_semantics_not_confirmed: call ConfirmAIDesignSemantics(object, design = list(batch = 'batch', ...)) before runni |
| 6 | plan_shape | plan | value_reference_unresolvable: action step_7 references ${step_6.output.reduction} (step step_6=run_integration exposes: class, n_features, n_cells, assays, redu |
| 7 | plan_refused | plan | invalid params for action s6: unknown parameter(s): note \| action s6 references undeclared dependency step(s): s3; add them to this step's depends_on \| prerequi |
| 8 | plan_shape | plan | success_criterion_unresolvable: c9 references unknown step step_10 |
| 9 | plan_refused | plan | value_reference_malformed: action step_7 uses ${step_6.output}; write ${step_6.output.<field>} instead (the syntax is ${step_id.output.field}) \| success_criteri |
| 10 | plan_refused | plan | action step_11 depends on unknown step(s): step_10 \| success_criterion_unresolvable: c9 references unknown step step_10 \| prerequisites not met for action step_ |
| 11 | plan_shape | plan | action is not registered: not_run; the `action` field must be an exact registry key, and a placeholder such as "none" is not a way to drop a step -- remove the  |
| 12 | plan_empty | plan | plan_empty: the plan has no actions. A plan with zero steps performs no analysis and cannot be executed; re-propose the plan, or if the goal genuinely needs no  |
| 13 | plan_refused | plan | prerequisites not met for action step_7: reduction is not available: ${step_6.output.reduction} |
| 14 | plan_knowledge | plan | prerequisites not met for action round2_step_1: contrast_invalid: ident.1/ident.2 must be levels of 'condition': ctrl, stim |

## 结论（自动生成，按 spec 要求的句式；不允许写成单一平滑百分比）

> **m12a-batch-5947e4b.jsonl**：12 次 measured 里 11 次完成（91.7%，95% 区间 [64.6%, 98.5%]）；失败集中在 plan_knowledge（1 次）；commit：5947e4b。

> **realdata_kang_ifnb_bench.jsonl**：42 次 measured 里 28 次完成（66.7%，95% 区间 [51.6%, 79.0%]）；失败集中在 plan_empty（2 次）、plan_shape（3 次）、plan_knowledge（1 次）、plan_refused（6 次）、execute_failed（2 次）；commit：4400568, fc3b7ff, a213a38, 60afa7c, 98aa395, 8d88a3b, b91ac54, ef8e4e6, 2c52dbb, 04c20aa, 80caa2e, d174367, 3e4ed58, c84a899, a8b1cc7, 2bf7b4f, 15fe726, 568a1da, 5947e4b。

往下一轮移交的东西：**失败形状计数表**。修掉其中一个形状，就重跑一遍，
看它是否真的从表里消失——那是唯一能证明修复生效的证据。

### 读这张表时必须同时说的话

- 区间宽度就是不确定度本身。12 次得到的 [64.6%, 98.5%] 与 42 次得到的
  [51.6%, 79.0%] 重叠，所以「新树更好」是**方向性证据，不是已证实的差异**；
  要分开两者需要更多次，或跑同一批 fixture 的对照。
- 跨 commit 的汇总（第二行）**不能**当作当前树的成功率：它里面有已修的缺陷。
  当前树的数字只有标注了单一 commit 的那一行算。
- `report_unverifiable` 为 0 **不等于**自检从未失败：M12-e 之前没有渲染步骤，
  这个桶还没有内容可装。
