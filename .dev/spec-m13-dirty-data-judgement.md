# Spec: M13 — 脏数据上的新判断，以及 D3 的机读实现

Status: draft（2026-10-21 起草）
Baseline: `3bb012b`（`make check` 0/0/0；AI 套件 FAIL 0 / PASS 5951）
上游: `.dev/spec-m12-plan-quality-and-reliability.md` §12 是本节的前身；它当时只固定了
方向，本节补细节。§12 的分歧点已在下面重新核对（见 §2.2）。

---

## 1. 这一轮要回答的唯一问题

> **给一份从没见过、且脏的真实数据，agent 会不会做出 M12 之前不会做出的判断？**

注意问法。不是"成功率会不会更高"——成功率高可能只是把失败救回来，产出完全一样。
**新判断才是能力。** 这一句决定了整个 spec 的验收方式：判据必须是"计划里出现了
M12 之前不可能出现的东西"，而不是"跑通了几次"。

M12 结束时的真实边界（不是谦辞，是可核对的）：

| 能说 | 不能说 |
|---|---|
| 同一 goal 12 次里 11 次完成，区间 [64.6%, 98.5%] | 面对没跑过的数据能想对 |
| 四维判据里 D1/D2/D4 能算 | D3（阶段递进）能算 |
| 六题（M1–M11 已知答案）6/6 | 没见过的局面 |
| `data_support` 让"数据还能做什么"可读 | 模型因此规划得更好（**从未在线验证过**） |

最后一行是最要紧的：M12 证明了 **`data_support` 字段是对的**，但**从未证明模型会用它**。
M12-b 报告 §6.4 自己写明了这一点。M13 必须补上——否则 M12 的收益只存在于账本里。

---

## 2. 范围

### 2.1 四个子件

| 编号 | 子件 | 回答什么 | 预估 |
|---|---|---|---|
| M13-a | 脏数据闭环 | §1 的主问题 | 1.5–2 天（含数据准备） |
| M13-b | D3 机读实现 | 「这轮因上轮说了什么而改变」能否变成数 | 1 天 |
| M13-c | `data_support` 的在线效用 | 模型是否真的读它、并因此改变计划 | 0.5 天（与 M13-a 共用采样） |
| M13-d | 验收物与移交 | 让 M13 的结论能被引用 | 0.5 天 |

### 2.2 与 §12 分歧点的重新核对（本节是上一节承诺的）

§12.4 建议「优先选现有 gate 会拒绝、而正确做法不是拒绝的脏形状，因为拒绝已经做得很好」。
**重新核对后，这个标准要改。** 理由如下：

- 「gate 会拒绝」的形状确实已经测得不错，但那批拒绝都是**数据缺失型**（没有 counts、
  没有 spatial），正确答案是「不规划」，**不产生任何新判断**；
- 真正对应「新判断」的是另一种：**数据没缺，但有一个正确做法是"换方法/换顺序/先清理"**。
  ambient RNA 就是这个形状（数据齐，但不清就没法解释）。

所以 §2.3 的选型标准改为 §12.5 的下位版本：**答案必须是"换个做法"，不能是"拒绝"**。

---

## 3. M13-a　脏数据闭环

### 3.1 数据选型

候选四种，逐一核对可得性后选定。判据：**答案必须是"换做法"，不能是"拒绝"**。

| 候选 | 难在哪 | 正确答案的形状 | 是否可得 | 选不选 |
|---|---|---|---|---|
| **ambient RNA 主导** | 空滴/低质量细胞的 RNA 混在真实群边缘，双细胞与环境 RNA 难分 | **先做 ambient 校正，再解释群组成**——不是"不解释"，也不是直接解释 | 需下载（见 §3.2） | ✅ **首选** |
| 多条件 × 深度混杂 | condition 与测序深度相关 | 把混杂写成限制——M2 已经做过一次 | Kang 已有 | ❌ 与 M2 重复 |
| 超稀疏（低 depth） | marker/velocity/milo 全部失真 | 换模型/换阈值——**但判据难以客观化** | 需构造 | ❌ 不好验收 |
| 肿瘤内异质性 | 亚群与拷贝数状态重叠 | 先 CNV 再谈亚群——M10 已经做过 | GSE144735 已有 | ❌ 与 M10 重复 |

**ambient RNA 胜出的理由**：(1) 正确答案是一个明确的**动作**（`run_decontx`），
使"新判断"可被判定——模型要么规划了它，要么没有；(2) 它是现有六题 Q1–Q6 里
**唯一没有覆盖**的主要质量问题（Q1 是 doublet，不是 ambient）；
(3) 包内 action 已有对应响应（`RunDecontX`），不需要新能力。

### 3.2 数据来源与落盘

**选择：GSE158055 COVID PBMC**（已在 `/home/wang/data/source/omics/data/covid_pbmc_gse158055/`，
含 cell annotation、counts、特征表、sdrf）。

需要做的只是把它变成 fixture，与 Kang 同一形态：
`R/fixtures_covid_pbmc.R`（新建），参照
[`R/fixtures_kang_ifnb.R`](../fixtures_kang_ifnb.R) 的既有模式：

- `covid_pbmc_paths()` / `covid_pbmc_download()`（文件已在则跳过）；
- `fixture_covid_pbmc()`：读 counts + annotation，转 SCE，**缓存成 rds**（同
  `kang_ifnb.rds`，48.9 MB 量级）；
- `fixture` 名与尺寸记进每条结果行，同既有约定。

**必须先跑一次 smoke**：`ncol` / `nrow` / anno 列 / 质量分布打印出来，
确认它真的是 ambient 主导（否则这一轮的前提不成立，§3.1 的判据失效）。
**如果 smoke 显示 ambient 不显著，换候选并重写 §3.1，不许硬跑。**

### 3.3 判定方式：`data_support` 前后对照

同一份数据、同一个 goal，跑两次计划：

| 组 | 上下文 | 预期 |
|---|---|---|
| **对照组** | 屏蔽 `analysis_story$data_support`（`options` 或上下文注入开关） | 可能直接规划 marker/注释，**不做 ambient 清理** |
| **实验组** | 正常注入 `data_support` | 应规划 ambient 清理，或在 rationale 里说明为什么不清 |

**注意这两组的做法**：不是"跑通没有"，而是**计划里有没有出现 ambient 清理这一步，
以及它出现在什么位置**。

为了让"新判断"可判定而不是主观印象，需要一条**确定性的前置判据**，从数据本身算出来，
不依赖模型：

```
ambient_pressure: 由数据计算的标量，例如
    (全局 zero-rate 与 top 表达基因 zero-rate 的偏离) × (总 UMI 的分布下尾比例)
```

这个判据只回答一件事：**这份数据是否确实存在 ambient 压力**。它由 R 确定性计算，
与被测模型无关。判据的实现与阈值要先在 §3.2 的 smoke 上标定，并把标定结果写进报告。

### 3.4 不做什么

- **不新增任何 action。** `RunDecontX` 已存在，缺的是模型知道该用它；
- **不评价 DecontX 的生物学正确性**（同 M5 不测注释准确率）；
- **不追求跑满整个闭环**。验收物是「计划 + 报告」，不是「一条完整的单细胞分析」。
  这一轮测判断，不测覆盖。

---

## 4. M13-b　D3 的机读实现（B18）

### 4.1 定义（沿用 M12-b §4.2 的 D3）

> 第二轮之后的计划，是否**真的读了上一轮的报告**并据此改变。

四条可判定判据（都从两份计划的 diff 得出，**不读模型内心**）：

| 编号 | 判据 | 计算方式 |
|---|---|---|
| D3.1 | **回应上轮告警** | 上轮报告的 `warnings`/`limitations` 逐条比对：本轮是否 (a) 有对应动作消除它，或 (b) 在 rationale 显式携带它并说明理由。**都不做 = 未回应** |
| D3.2 | **无进展检测** | 本轮的 action 集合与上轮的交集，且参数未变的，占全部 action 的比例。比例 > 0.5 → `no_progress` |
| D3.3 | **证据缺口闭合** | 上一轮报告 `evidence_gaps` 里每一条，本轮是否有动作对应 |
| D3.4 | **用户决策被沿用** | 上一轮记录的 `user_decision`，本轮是否真的使用了对应值（澄清通道的第二轮问题） |

### 4.2 实现要点

- 输入是**同一对象上相邻两轮**的 plan 与 report；计算放在
  `sclet_ai_round_progression()`（新内部函数，`R/ai-data-support.R` 或新文件）。
- 输出必须能被 `notes/` 里的记录直接引用，所以字段要平：
  `n_warnings_answered` / `n_warnings_ignored` / `progress_ratio` /
  `gaps_closed` / `decisions_used` / `progression_class`
  （`advanced` / `confirmatory` / `no_progress` / `regressed`）。
- **`regressed` 是必须存在的一类**：本轮撤销了上一轮已确立的结论（例如重新用一份
  已被标记为污染的 cluster 群做解释）。没有这一类，D3 就只能说好与一般，说不出坏。
- `progression_class` 进入 M12-e 的 `report$verification_status` 旁边，让**下一轮
  决定要不要继续时看得见上一轮算是什么**。

### 4.3 验收物

`test-ai-round-progression.R`：四类各至少一个合成 fixture（`advanced` /
`confirmatory` / `no_progress` / `regressed`），断言分类与计数正确；
外加一条"同一对计划两次调用结果一致"的确定性 pin。

### 4.4 一个必须写明的限制

D3 只能对**同一对象上真实存在的两轮**计算。它**不能**对单轮上下文生效，所以
`AIPlanAnalysis()` 在第一轮拿不到 D3——**这是设计而非缺陷**：第一轮没有前一轮。
但报错信息必须说清是「第一轮，无前序」而不是「D3 不可用」，否则读者会以为算错了。

---

## 5. M13-c　`data_support` 的在线效用

这是 M12 欠的那一笔（M12-b 报告 §6.4）。

### 5.1 方法

复用 M12-a 的 12 次采样脚本，但把**同一个 goal 各跑一次（对照）与正常（实验）**：

- **对照组**：把 `data_support` 从上下文里剥掉（不改代码，用一个显式的
  `options(sclet.ai.omit_data_support)` 或上下文参数），跑 N 次（N 与 M12-a 相同，12）；
- **实验组**：正常跑，用 M13-a 的脏数据。

统计的不是成功率，是 **`unexplored` 里的建议有多少被真正规划**：

```
adoption_rate = |plan_actions ∩ unexplored_actions| / |unexplored_actions|
```

**这个数是 M12 收益的直接度量**：如果 `data_support` 是对的但没人用，adoption ≈ 0，
那 M12 的结论要改写（"字段对，但模型不用"），而不是继续声称收益不明。

### 5.2 一个诚实的陷阱

「adoption 低」可能有两种原因，必须分开，不许混成一个数：

1. **模型没读到**（注入失败、topic 选择漏了）→ 可修，修的是接线；
2. **模型读到了但认为不该做**（`unexplored` 建议了，但目标不需要）→ **不是缺陷**，
   `unexplored` 是数据推断，目标是人的意图。

所以每条 trial 要同时记录「是否被注入」（`sclet_skills_enabled` 已有此字段，
参照它加一个 `data_support_injected`），否则 adoption 低的原因无法归因。

### 5.3 验收物

`notes/m13-data-support-adoption.md`：两次采样的 adoption 对比 + 归因分栏
（没读到 / 读到了但没采用 / 读到了且采用了）。

---

## 6. 执行顺序

**b → a+d（并行）→ c → 验收**

1. **M13-b 先做**（1 天）：它是纯 R 确定性代码 + 单测，不依赖在线调用，
   且 M13-a/c 的验收要用它的输出；
2. **M13-a（数据 + 闭环）与 M13-d 并行**：M13-a 要先完成 §3.2 的 smoke 才有意义，
   smoke 失败就得换候选——所以这一子件的**第一步是验证前提，不是跑分析**；
3. **M13-c 最后**（0.5 天机时，约 36 次调用：对照组 12 + 实验组 12 + 已有），
   因为它要同时用上 M13-a 的 fixture 和 M12-a 的对照方法；
4. **M13-d 收口**：写验收物，并更新 §2.2 那张表（能说/不能说）。

### 6.1 明确的顺序理由

M13-b 排第一，因为**它是唯一能独立于模型质量完成的部分**。如果 M13-a 因数据问题
延期或换候选，M13-b 不受影响；反过来先做 M13-a，一旦换数据，b 的验收数据也要换。

---

## 7. 不为 M13 做的事（写给下一个里程碑）

- **注释准确率**：继续不测（M5 决定）；
- **spatial / multimodal 转正**：需要对应类型的数据，是 M14 量级；
- **多 provider 横向对比**：会和 provider 变更混在一起；
- **新一轮 offline 题库（Q7+）**：本轮的主判断在预见过的之外，添加越多已答案的题库，
  只会稀释"新判断"这一指标——除非某个新形状被真的发现了，那时自然要补。

---

## 8. 风险与预设处置

| 风险 | 处置 |
|---|---|
| GSE158055 的 ambient 压力不显著 | §3.2 smoke 显式检查。不显著就**换候选**，按 §2.2 的新标准重选，不许硬跑 |
| ambient_pressure 判据标定不出来 | 退化判定：只要有 `run_decontx` 出现，记为采纳；没有则记为未采纳。**并说明阈值未标定**，不假装客观 |
| 对照组与实验组 success rate 差异被误读 | 本轮的指标不是成功率。报告里成功率只作为背景，**adoption 与 D3 才是主判据** |
| adoption 低 | §5.2 已规定必须归因到「没读到」还是「读到了不用」。归因不明 = 本轮失败 |
| M13-b 的 `regressed` 判不出来 | fixture 补一个显式回归用例；如果真实数据里从没出现，就如实写"未观测到"，不编 |
| 成本 | M13-a 约 12 次 × 6 min ≈ 72 min；M13-a smoke + 数据准备约 1 h；M13-b 无调用。总计一次付费采样 |

---

## 9. 一句话

> M12 证明"可靠"可以测量；M13 要证明"会分析"可以测量——**用一份没跑过的脏数据，
> 看计划里是否出现了 M12 之前不可能出现的判断**。同时把 M12 欠的一笔补上：
> `data_support` 字段是对的，但模型到底用不用，本轮给一个有归因的数。

---

## 10. M13-b 实施记录（2026-10-21）

### 10.1 交付物

`sclet_ai_round_progression()`（`R/ai-data-support.R`，新）+ `sclet_ai_previous_round_progression()`。
D3 现在**有数了**，四个类：`advanced` / `confirmatory` / `no_progress` / `regressed`，
外加 `first_round`（首轮，不是缺陷）。

`report$round_progression` 已接进 `RunAIAnalysis()` 的返回对象。
pin：`tests/testthat/test-ai-round-progression.R`，**43 断言**。

### 10.2 实现中遇到的三个真实问题（都不是预设的，写下来以免重复）

**问题一：`return()` 在 `tryCatch` 里穿透外层函数——第三次。**

这个包里的第三次。`sclet_ai_previous_round_progression()` 里我在 tryCatch 表达式内写了
`return(NULL)`，结果从**外层函数**返回裸 NULL，D3 全部字段消失。
（第一次是 M12-b 的 `has_replicates`，第二次是同一轮的其它 early exit。）
修法同样是赋值 + 值传递，不用 `return()`。

**这条要升级成代码规范**：一个 helper 里出现 `tryCatch({ ... return(...) ... })`
就必须改。它已经在三个不同函数里以三种不同伪装出现过两次，还会再来。

**问题二：上一轮不是一个记录，而是两个。**

D3 需要"上一轮的计划"和"上一轮的报告"。**没有任何单个记录同时含两者**：

- 计划（跑了什么）在最新 `ai_execution` 记录的 `params$steps` 里；
- 报告（于是说了什么）在最新 `ai_report` 记录里，且其 `warnings` 在
  **`artifacts$warnings`**（list of character），`summary$cited_user_decisions` 才是有引用的用户决策。

我最初的实现去 `ledger$state_records$ai_ai_report` 找报告、去 `summary$warnings` 找告警，
**在两个真实对象上都得到 `not_comparable`**。而这个失败读起来跟"D3 从未被计算过"一模一样——
如果那时我只跑单测得出的结论会是"D3 实现有问题"，而真相是**我找错了地方**。
所以那条 pin 用了真实数据 probe（跳过保护），这是它唯一的验证途径。

**问题三：两处 NULL/定长崩溃。**
`nzchar(NULL)` 是运行时错误；`vapply` 要求 FUN 定长返回，而 topic helper 返回 NULL。
两个都在真实数据上崩过（合成 fixture 恰好都绕开了）。都已修，各有一条 pin。

### 10.3 一个值得单独说的结果

在 m7 对象上实测，D3 给出 **`regressed`**，而且它是正确的：

上一轮（M7 报告）的 warning 里有
"cluster 12 的 Monaco 标签（单核）来自被双细胞主导的群，报告组成时应剔除或标注"，
而当前计划重新做 `run_annotation` 与 `run_de_test`，**没有先做 decontamination**。
这正是 `regressed` 这一类要抓的东西——它重新解释一个已知被污染的群，
比什么都不做更贵。

**这个判断是引擎算出来的，不是人写的**；而它在 M7 自己的报告里已经被人工记录过。
两条独立路径得到同一结论，这比多跑十次计划更能说明 D3 判据抓对了形状。

### 10.4 与 §4 的偏差

| §4 说了什么 | 实际 |
|---|---|
| D3.1–D3.4 四条判据 | ✅ 四条都实现，`no_progress` 用参数相等判定（同 action 同 params 才算重复） |
| `progression_class` 进入 report | ✅ `report$round_progression` |
| 「第一轮要报 first_round 而非 unavailable」 | ✅ §4.4 的要求，实现为 `first_round` |
| 首轮拿不到 D3 是设计而非缺陷 | ✅ 已实现并在测试中固定 |

---

## 11. M13-a 实施记录（2026-10-21）

### 11.1 数据与前提（都验证过，不是假设）

- **GSE158055 单样本 38**：27,943 genes × 8,129 cells，15,420,687 nonzeros。
- 全队列 146 万细胞 / 23 亿 triplet 被**明确拒绝**：规模使 decontX 不可跑，
  且混样本会让判据测到批次而非 ambient。理由写进 fixture 注释。
- **提取只扫了文件的 12%**：barcode 按样本连续分组 + mtx 按列排序，扫到目标样本末列即停
  （2.79 亿 / 23 亿）。提取器有**连续性断言**：样本列不连续就大声失败，不返回残缺矩阵。
- **前提成立**：`ambient_pressure` = **high**（low-UMI tail 0.0335 ≤ 0.05；
  top200 检出率 0.9139 ≥ 0.90，两项都触发）。

### 11.2 一个推翻了本文档前提的发现

§3.1 原写「包内 action 已有对应响应（`RunDecontX`），不需要新能力」。**这是错的**：
`RunDecontX()` 导出但**没有 AIAction**。正确答案在当时不是可规划的一步。

按用户决定 **A + B 并做**：建 `run_decontx`（A），同时保留「如实承认做不到」这一支（B）。

### 11.3 A：`run_decontx`（33 → 34 个 action）

`R/ai-registry-preprocess.R`。要点：整数 counts 必需（非整数报 `counts_assay_not_count`）、
`celda` 必需、写 `decontXcounts` 与 `decontX_contamination`、**不删细胞、不切 active assay**、
证据节点只含分布且标 `contamination_call_is_model_estimate = TRUE`。
pin：`test-ai-decontx.R`，42 断言。

**实现中撞到四个问题，都修了：**

1. `source` 必须是真实记录 id（`"decontX"`），不是证据节点名——每次
   `source must identify a completed analysis`。
2. 证据值里不能有字符串：`cleaned_assay = "decontXcounts"` 被隐私层拒为无界文本。
   规则是对的，结构事实属于 methods 区不属于 finding。
3. `type = "preprocess"` 不等于做过 decontX（归一化也是 preprocess），判断必须按 id。
4. **decontX 是确定性的**：两次同输入污染向量完全相同（max abs diff = 0），
   所以重跑 no-op 是**正确**行为。这条尤其值得记：直觉会假设随机估计器会变，
   于是写出一条钉在错误前提上的测试。输入变了才应冲突，那一支单独 pin。

### 11.4 真实运行

```
elapsed 68.1 s | mean 0.0487 | median 0.0202 | q95 0.1735 | max 0.9915
evidence: ev:decontx_decontx_covid38
```

`data_support` 在该对象上只给一条建议，且是 `run_decontx`，`why` 带 screen 实测数字。

### 11.5 新增的确定性判据

`sclet_ai_ambient_pressure()`：**两个 screen 分开报、不相乘**（相乘会得到一个看不出
哪半动了的数），并把自己的局限写在返回值里（高检出率也是管家基因的样子）。
阈值在当前 fixture 上标定，**不是文献值**，代码里明写。

### 11.6 与 §3 的偏差

| §3 说 | 实际 |
|---|---|
| 判据只回答「有没有 ambient 压力」 | ✅ 两个 screen + label，且自带 caveat |
| 「不新增任何 action」 | ❌ **破例新增 `run_decontx`**，经用户明确同意。原因见 §11.2 |
| 验收物是「计划 + 报告」 | **部分**：本轮给的是前提验证 + 上下文 + 动作 + 真实运行；**「模型会不会用」留给 M13-c** |

报告：`notes/m13-ambient-report.md`。

---

## 12. M13-c 实施记录（2026-10-21）

### 12.1 方法

`sclet.ai.omit_data_support` 开关（新增，命名上就表明「关掉它上下文更差」，只服务对照臂）。
对照组剥掉 `data_support`，实验组正常，各 6 次，只规划不执行。
goal **刻意不提污染**——否则测的是「听不听话」而不是「会不会用数据推断的上下文」。

**必须先验的一项**：注入是否正确。12/12 全部符合预期（对照 FALSE、实验 TRUE）。
没有这一项，「adoption = 0」可能只是开关没生效。

### 12.2 结果

| 组 | 规划 `run_decontx` | 比率 | Wilson 95% |
|---|---|---|---|
| 对照 | 0 / 6 | 0% | [0%, 39%] |
| 实验 | 2 / 6 | 33.3% | [9.7%, 70%] |

**Fisher p = 0.455，不显著。** 诚实的读法：方向一致，但 12 次样本无法与偶然区分。

**对照组基线为 0 这一点比 33% 更有说服力**：`run_decontx` 已是注册动作，模型本来
可以自己想到它，6 次一次都没有——说明它不是模型的默认习惯。

### 12.3 真正的发现：4 次「拒绝」不是拒绝

未采纳的 4 次，rationale 里**一次都没提到 ambient**。它们不是 spec §5.2 里
「读到了但有理由地拒绝」（那不算缺陷），而是**建议在上下文里、却没进入推理**。

这与 M12 修的 topic 选择同类：信息有了，但没被捡起来。记为 **B20**。

另一条：t6 规划了 `run_decontx` 却没有任何解释——**做对了也要能说出来为什么**，
记为 **B21**。

### 12.4 对 §5 的偏差

| §5 说 | 实际 |
|---|---|
| 指标是 adoption 而非成功率 | ✅ |
| adoption 低必须归因「没读到」/「读到了不用」 | ✅ 并发现了第三类：**读到了但没进入推理**，这是最有用的那条 |
| 每 trial 记 `data_support_injected` | ✅ 12/12，且加了 `injection_matches_arm` 逐条校验 |
| N 与 M12-a 相同（12） | ❌ **每臂 6**（共 12）。理由：单次 40 s，两次全量更贵，且 6+6 已足以看出「是否有区分度」。代价是**统计上不显著**，报告里已明写 |

报告：`notes/m13-data-support-adoption.md`。
