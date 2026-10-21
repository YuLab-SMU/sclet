# M12-b 计划质量报告：`data_support` 上线前后的对照

> 状态：M12-b 完成（2026-10-21）
> 验收物形态：**规则（四维判据）→ 实现（`data_support`）→ 六个真实局面的实测结果 → 与 M1–M11 已知答案的对照**
> 脚本：`.dev/m12-question-bank.R`（六题，可单独跑）；pin：`tests/testthat/test-ai-data-support.R`（63 断言）

## 1. 这份报告要回答什么

M12-b 的主张是：**计划质量可以从账本读出来，而不靠观感。** 这份报告就是那句主张的证据。
它不测模型，测的是**上下文有没有把「正确答案」变成模型可读的东西**。

判据分两半，严格分开：

| 半边 | 内容 | 如何验证 |
|---|---|---|
| 机读 | 四维判据 D1–D4 里可计算的项（数据形状、gate 理由、已有分析集合、`unexplored` 建议） | 脚本 + pin 测试 |
| 人工 | 「这个建议在科学上对不对」 | 本节 §4 逐条对照 M1–M11 的既定结论 |

**人工判定不进统计**，这是 spec §4.2 定的规矩：一旦判据变成「看起来合理」，这轮就退化成写观感。

## 2. 实现了什么

`analysis_story$data_support`（新字段，schema 1.0；`analysis_story` 由 1.0 → 1.1）。
纯 R 确定性计算，不含模型调用、不含逐细胞值、不含受保护标识符。

| 字段 | 内容 | 对应判据 |
|---|---|---|
| `shapes` | 数据形状事实：counts 是否存在**且为整数**、spliced/unspliced、多样本、每条件重复、spatial、多模态、候选设计列 | D2 前提 |
| `supports` | 每个问题类别「数据能不能回答」，理由用**既有 gate 名**，附 `unblocked_by` | D2 |
| `blocked` | 同一信息的 gate 优先视图 | D2 |
| `unconfirmed_semantics` | 需要用户回答才能规划的选择（设计角色、注释参考） | D2 / 澄清通道 |
| `unexplored` | **数据能做、账本里没有**的分析：`what` + `action` + `why` | D4 |

`unexplored` 是唯一没有前身的一项：`evidence_gaps` 说「某分析缺证据」，
`unexplored` 说「**这份数据从来没被分析过**」。

## 3. 六个真实局面的实测结果

数据来自 M1–M11 的既有对象（`/tmp/sclet-agentrun/`），非重建。

| 对象 | 细胞数 | 多样本 | 重复 | 整数 counts | spliced | `unexplored` | 建议内容 |
|---|---|---|---|---|---|---|---|
| m2-final（Kang） | 43,095 | ✅ | ✅ | ✅ | — | 4 | milo, trajectory, enrichment, annotation |
| m5-report（PBMC 注释） | 11,769 | — | — | ✅ | — | 1 | trajectory |
| m7-report（PBMC 双细胞） | 11,769 | — | — | ✅ | — | 1 | trajectory |
| m10-report（CRC CNV） | 1,400 | ✅ | — | ✅ | — | 1 | trajectory |
| m11-report（pancreas 速度-命运） | 3,696 | — | — | ✅ | ✅ | 2 | markers, annotation |
| m1-report（PBMC 主干） | 11,769 | — | — | ✅ | — | 2 | trajectory, annotation |

**六个对象的建议逐条核对，无一条是「已做过的事」。** 最有说服力的是 m11：
它做过速度、轨迹、命运，`data_support` 正确地说「velocity 和 fate 都别重跑」，
同时指出两件从未做过的事——**marker 检验和参考注释**。这正好是「探索数据价值」要的输出。

## 4. 与 M1–M11 已知答案的对照（回归题库）

六题，每题都有一份**已记录的正确答案**，判据从那份记录来，不是事后编的。

| # | 局面 | 已知正确答案（出处） | 判据 | 结果 |
|---|---|---|---|---|
| Q1 | PBMC cluster 12 是 98.3% 双细胞（M7） | 排除双细胞解释，且不重复已做的 doublet | 不重跑 doublet **且** 提出尚未做的轨迹 | PASS |
| Q2 | velocity 需要 counts，本仓库导出只有 logcounts（M9） | 拒绝 velocity，理由 `counts_assay_missing` | 同上 + 前提自证（对象确实缺整数 counts） | PASS |
| Q3 | Kang 的 stim 只在 batch2（M2） | 把设计混杂当限制，识别多样本与 batch 列 | 识别多样本 + batch 列 + 重复 + 未确认语义 | PASS |
| Q4 | CNV 签名阈值取中位数 → 恰好一半（M10） | 不重跑 CNV；**0 条建议必须有理由** | 不重跑 + 空清单须因 `reduction_missing` 成立 | PASS |
| Q5 | cluster 9 主导标签 0.22 → 混合，不给类型（M5） | 注释已完成，不重复 | 不重跑 annotation **且** 清单非空 | PASS |
| Q6 | 终态占比 vs 去向分布（M11） | 不重跑 velocity/fate | shapes 正确 + 不重跑 + 清单非空 | PASS |

**6/6 通过。**

## 5. 题库抓出的两个真缺陷（这是本轮最有价值的部分）

题库第一版有六题全 PASS 的结果，但**那是假的**。加强判据后暴露两个缺陷，都已修并 pin：

### 缺陷一：空清单被误判为正确（第一版判据太弱）

第一版的 Q1/Q4/Q5 只断言「不建议重跑已有分析」。一个**永远返回空 `unexplored` 的实现**
能轻松通过——而这正是「空结果被当成正确结果」那类缺陷，恰好是可靠性定义要防的事。
改成同时断言**清单非空且指向真实缺口**（Q4 例外：允许为空，但必须证明空得有理）。

### 缺陷二：把记录类型名当成分析存在（真代码缺陷）

**m11 的 `unexplored` 当时是 0 条**，判 FAIL。核查后确认这是代码错，不是数据特殊性：

- pancreas 对象有 `preprocess`、`reduction`、`velocity`、`trajectory`、`cellrank`，
  **但没有 `clustering` 记录**——那次运行从嵌入直接走到轨迹，中间没有独立的聚类步骤记录；
- 我的代码用 `has_analysis_type("clustering")` 判断「有没有分组」，于是判成
  `cluster_identity_missing`，并让所有依赖分组的建议全部失效。

**这正是 AGENTS.md 里「Message–artifact disagreement is a defect class」的同构**：
记录类型名是 artifact 的标签，不是分析是否发生的证据。修法是新增
`sclet_ai_data_support_grouping_evidence()`：**任何以分组为前提的分析都是分组的证据**
（detest、enrichment、annotation、trajectory、velocity、cellrank、milo…），
外加对象自带的 `sclet_ident` / `label`。修复后 m11 从 0 条变成 2 条，
而且这两条正是它真正缺的。

### 附带：Q4 的「0 条」是正确答案，不是缺陷

m10 是 inferCNV 的原始对象，**没有降维也没有聚类**，所以轨迹那条被
`reduction_missing` 挡住是**对的**。修复缺陷二之后它变成 1 条建议（轨迹仍被挡，
但别的缺口浮现）。这条单独立出来，是因为**「0 条」需要被证明，而不是被假设**。

## 6. 没做到 / 要如实说的

1. **D3（阶段递进）没有机读实现。** 它要求「本轮计划是否回应了上轮的 warnings/limitations
   并据此改变」，需要**两轮计划**才能计算，而 M12-b 目前只处理单轮上下文。
   本轮只把它写进了 `sclet-analysis` 技能（模型可读），**没有变成可算的数**。
   这是与 spec §4.2 的偏差，记录在此，不掩盖。
2. **`unexplored` 的完整性只覆盖 9 类**（integration、milo、velocity、trajectory、
   markers、enrichment、annotation、spatial 等）。它不试图穷举——一个「什么都建议」的清单
   会让模型失去判断力。**漏掉的建议是保守的沉默，不是错误**；但要说清这是设计选择。
3. **六题全部是「见过的局面」。** 这正是 M13 要补的：题库检验「面对见过的局面能不能想对」，
   M13 检验「面对没见过的局面能不能想对」。**不要把 6/6 读成能力证明。**
4. **`data_support` 尚未在真实模型调用中验证。** 它进入了 `AIPlanAnalysis()` 的上下文，
   但本轮的证据是「字段算得对」，不是「模型因此规划得更好」。后者需要 M12-a 的在线重复测量。

## 7. 对 §4.4 验收物四条的回答

| 要求 | 状态 |
|---|---|
| 同一数据同一目标，人工计划与 agent 计划并列 | **部分**：本轮给的是判据 + 六题对照，未做「人工 vs agent」双栏（需要在线调用，属 M12-a） |
| 四维逐条判定，机读项进表 | **D1/D2/D4 有机读实现；D3 只有技能层文字**（见 §6.1） |
| 差异说明逐条对应 step | 未做（同第 1 条，需在线计划） |
| `data_support` 上线前后 D2/D4 是否改变 | **有改变**：m11 从 0 条建议 → 2 条真实缺口；m10 从 0 → 1；六个对象全部给出可核对的具体建议 |

**结论一句话**：`data_support` 让「数据还能做什么」从模型的猜测变成可核对的字段，
六个真实对象上逐条正确、无一条建议已做过的事；但**它的效果证明只到「字段对」，
「模型因此变得会规划」要靠 M12-a 的在线测量来回答。**
