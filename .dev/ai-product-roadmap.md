# sclet AI 产品与架构路线图

> 状态：路线收敛版（2026-10-06）  
> 基线：`devel`，已完成 context / privacy / planning skeleton / controlled execution / evidence / claim ceiling / integration / annotation / rare-cell / trajectory 第一批。

## 0. 这份路线图解决什么问题

sclet AI 的目标不是不断增加可以被 AI 调用的分析函数，而是让 AI 能够可靠地完成一条科研工作闭环：

```text
理解当前分析过程
    -> 识别缺失语义与风险
    -> 提出可验证的分析计划
    -> 在用户确认后执行
    -> 依据 evidence 解读结果
    -> 形成下一轮计划
```

后续工作必须优先服务这条闭环。新增领域 action 只有在不破坏闭环、且能复用统一 contract 时才进入路线。

## 1. 当前判断

### 1.1 已经稳固的基础

- **Bounded context**：`GetAnalysisLedger()`、`GetAIProfile()`、`sclet_ai_context()` 提供不含完整表达矩阵、原始 barcode、患者信息的分析状态视图。
- **受控规划与执行**：`AIPlanAnalysis()`、`ValidateAIPlan()`、`ExecuteAIPlan()` 支持 action allowlist、参数 schema、依赖、fingerprint、dry-run 和人工确认。
- **证据与主张约束**：`RecordAIEvidence()`、`ValidateAIEvidenceRefs()`、`sclet_ai_claim_ceiling()`、`RecordAIResult()` 提供 evidence refs、独立性和 claim ceiling。
- **隐私与设计确认**：provider-side privacy gate 和 `ConfirmAIDesignSemantics()` 已成为执行前边界。
- **已有领域 adapter**：integration、marker/DE/annotation、rare-cell/doublet、trajectory 第一批均已有 action 或 diagnostics。

### 1.2 当前的主要缺口

1. **分析历史仍主要是状态快照**：AI 能看到做过什么，但还不能稳定地讲清楚分析阶段、输入输出、关键决策、未解决冲突和证据缺口。
2. **计划 contract 偏 action-centric**：已有 action 和依赖，但假设、候选路线、成功标准、停止条件和风险还没有成为强制结构。
3. **结构校验强于科学合理性校验**：`ValidateAIPlan()` 能阻止非法 action，但对数据规模、科学问题和成功标准的适配性约束仍不足。
4. **解读没有成为执行闭环的一部分**：`AIExplainAnalysis()` 等入口存在，但执行完成后不会自动生成统一的 evidence-linked interpretation report。
5. **领域 action 扩展速度可能超过核心闭环建设速度**：velocity、CellRank、spatial、multimodal 继续扩展会扩大维护面，却不自动解决上述问题。

## 2. 产品原则

### 2.1 三层分离

所有 AI 能力必须明确属于以下一层：

1. **Facts / Context**：R 确定性生成的事实、状态、诊断和 provenance。
2. **Decisions / Plan**：用户确认的实验设计和 AI 提出的、待确认的分析计划。
3. **Interpretation / Report**：引用 evidence 的 AI 解读、主张强度、不确定性和替代解释。

AI 文字不能伪装成 deterministic fact；用户决定不能被 AI 重新描述成 AI 自己的判断。

### 2.2 统一用户流程

```text
Understand
  GetAnalysisLedger / GetAIProfile / deterministic diagnostics
      |
Clarify
  结构化问题 -> 用户回答 -> user_decision / design confirmation
      |
Plan
  目标、假设、候选路线、成功标准、停止条件、风险
      |
Validate
  schema / dependency / fingerprint / prerequisite / scientific sanity
      |
Preview
  dry-run，说明将执行什么、为什么执行、如何判断成功
      |
Execute
  人工确认 -> ExecuteAIPlan -> analysis state + evidence
      |
Interpret
  evidence-linked report：结论、不确定性、替代解释、限制、下一步
      |
Continue
  从 report 生成下一轮候选 plan
```

### 2.3 不把 action 数量当作进度

后续完成度以以下指标衡量，而不是以新增 action 数量衡量：

- AI 能否准确描述当前分析状态；
- 计划是否包含可验证的成功标准和停止条件；
- 执行前是否能发现缺失输入和不合理假设；
- 结果解释是否引用正确 evidence，并保留不确定性；
- 下一轮计划能否从上一轮报告自然产生；
- 用户是否能在不阅读内部 ledger 的情况下理解过程。

## 3. 统一 contract

### 3.1 Analysis Context

现有 context 的基础上，逐步增加 bounded `analysis_story`：

```text
schema_version
dataset
current_stage
active_view
completed_steps
analysis_timeline
inputs_outputs
user_decisions
assumptions
evidence_gaps
contradictions
blocked_actions
quality_checks
capabilities
privacy
fingerprint
```

要求：

- 由 R 确定性构建；
- 不包含完整 assay、barcode、患者信息或逐细胞原始值；
- 明确区分 observed state、candidate interpretation 和 user decision；
- 能指出“已经做了什么”“还缺什么”“哪些结果不能直接比较”。

### 3.2 Analysis Plan

在现有 `sclet_ai_plan` 上逐步补齐：

```text
goal
scientific_question
assumptions
design_semantics
candidate_routes
selected_route
actions
dependencies
expected_outputs
success_criteria
stop_conditions
risks
human_confirmations
```

要求：

- 每个 action 都有 expected output；
- 每个 success criterion 都能映射到可观测 metric、state 或 evidence；
- 变更对象的 action 前必须有必要的诊断或确认；
- 不能把用户尚未确认的语义写成已确认事实；
- 不能有“执行后看情况”这类不可验证的成功标准。

### 3.3 Analysis Report

新增统一 report contract，具体 API 名称可在 P0 决定，候选名称为 `AIInterpretAnalysis()` 或增强 `AIExplainAnalysis()`：

```text
schema_version
analysis_refs
summary
findings:
  - statement
    evidence_refs
    claim_level
    uncertainty
    alternative_explanations
    limitations
success_criteria_status
warnings
next_steps
proposed_followup_plan
metadata:
  context_fingerprint
  read_only
  actions_executed
```

要求：

- 每个事实性 finding 都要有 evidence ref；
- claim level 不得超过 evidence ceiling；
- 必须区分 measured / observed / associated / consistent_with / hypothesis；
- AI report 本身不是 deterministic evidence；
- 不记录模型隐藏思维过程，只记录输入 fingerprint、引用、结构化结论和用户决策。

## 4. 分阶段路线

### P0：路线冻结与 contract 对齐

**目标**：停止横向扩展，固定产品语言和接口职责。

**工作项**：

- 以本路线图为统一 roadmap；
- 对 `AIInvestigate()`、`AIPlanAnalysis()`、`ExecuteAIPlan()`、`AIExplainAnalysis()`、`RunAIAnalysis()` 做职责映射；
- 明确 context / plan / report 三种结构的边界；
- 为已有 integration、annotation、rare-cell、trajectory action 标注“输入、输出、证据、解释入口”；
- 不新增 velocity、spatial、multimodal action。

**验收**：

- 每个公开 AI API 只有一个主要职责；
- 文档中不再出现相互矛盾的“已实现/未实现”状态；
- 一条端到端 user journey 可以从 context 追踪到 report；
- 现有 focused AI tests 和 `make check` 保持通过。

### P1：补强 Understand / Analysis Story

**目标**：让 AI 能解释“这个对象经历了什么”。

**工作项**：

- 在 `GetAnalysisLedger()` 或独立 builder 中增加 bounded `analysis_story`；
- 记录分析阶段、时间顺序、输入/输出引用和 active state 变化；
- 汇总用户决策、设计确认、blocked actions 和 evidence gaps；
- 让 `AIInvestigate()` 按问题目标调用对应 deterministic diagnostics，而不是固定只提供 integration readiness；
- 对分析记录冲突、重复记录和缺失 provenance 给出明确 warning。

**验收**：

- 给定一个经过 Normalize -> PCA -> neighbors -> clusters -> annotation 的对象，AI context 能按顺序说明已完成步骤；
- 能区分“已测量结果”“用户确认语义”“AI 候选解释”；
- 能指出下一步的缺失输入，而不是泛化地建议“继续分析”；
- context 不包含逐细胞原始值和受保护标识符。

### P2：补强 Plan / Scientific Sanity

**目标**：让 AI 的计划不仅结构合法，而且可验证、可解释。

**工作项**：

- 扩展 plan schema，加入 assumptions、candidate routes、success criteria、stop conditions 和 risks；
- 在 `ValidateAIPlan()` 中加入确定性的 scientific sanity checks；
- 增加数据规模、输入 modality、已有 state 和目标 action 的兼容性检查；
- 要求先诊断、后有副作用执行；
- 将 clarification 结果直接转为 plan 的 human confirmations；
- 在 dry-run 中显示成功标准、停止条件和用户需要确认的项目。

**验收**：

- 结构合法但缺 success criteria 的 plan 被拒绝或明确降级为 draft；
- 缺少必要 assay/reduction/metadata 的 plan 被阻止；
- 未确认 batch/reference/root 的 plan 不能执行；
- dry-run 能展示完整 action、输入、输出、风险和验证方式；
- plan 执行失败后能说明是 prerequisite、action error 还是 success criterion 未满足。

### P3：补强 Interpret / Evidence-linked Report

**目标**：让执行结果自动进入可解释、可追溯的报告闭环。

**工作项**：

- 统一 `AIExplainAnalysis()` 与 `AIInvestigate()` 的 report 输出 contract；
- 增加 analysis result -> evidence -> finding 的引用校验；
- 把 success criteria 的实际状态纳入 report；
- 强制输出 uncertainty、alternative explanations 和 limitations；
- 让 `RunAIAnalysis()` 在执行成功后可选生成 report；
- 把 report 的 context fingerprint、analysis refs 和 user decisions 写回 ledger，但不写入隐藏 reasoning。

**验收**：

- 每个 finding 都能追溯到 analysis/evidence；
- claim ceiling 违规的 report 不能被记录；
- 零结果、缺失指标和冲突结果不会被润色成阳性结论；
- report 能明确回答“支持什么、不能支持什么、下一步做什么”；
- 从 report 生成的 follow-up plan 能通过同一套 ValidateAIPlan()。

### P4：稳定已有领域 adapter

**目标**：在扩展新领域前，证明现有 adapter 能完成真实任务闭环。

**优先验证对象**：

1. integration route comparison；
2. marker/DE/annotation；
3. rare-cell/doublet；
4. trajectory 第一批。

**工作重点**：

- 使用小型真实或高质量合成数据；
- 验证 context -> plan -> execute -> report 全链路；
- 检查 evidence refs 和 claim levels；
- 检查失败、部分成功、无结果和冲突指标；
- 收集用户真正需要的解释字段。

这一阶段原则上不新增领域 action，只完善已有 action 的输入、输出和解释接入。

### P5：受控扩展动态、空间和多模态

只有 P1–P4 完成后才进入。

建议顺序：

```text
velocity readiness
    -> velocity execution
    -> trajectory / velocity interpretation
    -> CellRank / fate
    -> spatial
    -> multimodal
```

每个新领域必须同时交付：

- readiness diagnostics；
- action adapter；
- prerequisite 和 design confirmation；
- success criteria；
- evidence summary；
- interpretation report 接入；
- 失败/不确定性/替代解释测试。

不能只增加一个 `run_*` action 就算完成。

## 5. 明确暂停的方向

在 P1–P4 完成前，暂停：

- `run_velocity` 的完整执行扩展；
- CellRank / fate action；
- spatial action catalog；
- multimodal action catalog；
- 大规模 pathway / communication action；
- evidence graph 的进一步复杂化；
- AI 文件的大规模移动或重命名；
- 仅为了增加测试数量而增加边缘 action。

以下内容不是暂停，而是必须持续守住的基础边界：

- privacy by default；
- dry-run + explicit confirmation；
- no arbitrary R/Python execution；
- no automatic cell deletion/filtering；
- no causal claim from observational evidence；
- no unconfirmed reference/batch/root semantics。

## 6. 第一条应真正跑通的示范闭环

建议用一个已有能力较完整的任务作为主线，例如 integration route comparison：

```text
用户问题：是否存在 batch effect？
    |
AI context：当前 assay、PCA、sample/batch 候选、已有分析、QC、缺失语义
    |
澄清：确认 batch、condition 和 preserve variables
    |
Plan：候选 fastMNN / Harmony / scVI、成功标准、风险、停止条件
    |
Validate + dry-run
    |
用户确认后执行
    |
Compare metrics + evidence
    |
Interpret report：支持什么、trade-off、未解决问题
    |
Follow-up plan：是否需要 marker/DE 或重新设计 integration
```

这个示范闭环比新增一个 velocity action 更能验证产品是否真的有用。

## 7. 路线决策规则

以后每提出一个新功能，先回答：

1. 它属于 Understand、Plan、Execute 还是 Interpret？
2. 它是否复用 context/plan/report contract？
3. 它是否能形成 evidence-linked 的可解释输出？
4. 它是否有明确的 success criteria 和失败处理？
5. 它是否能在不增加新的安全边界的情况下实现？
6. 如果不做这个功能，当前用户闭环是否真的被阻塞？

如果一个功能只是增加一个新的分析入口，但不能改善上述闭环，应暂缓。

## 8. 当前下一步

下一轮不直接开发 velocity 或 spatial。优先执行：

1. 将 `analysis_story` 和 context/plan/report contract 落到 spec；
2. 选择 integration route comparison 作为第一条完整示范闭环；
3. 设计 P1 的最小字段和测试 fixture；
4. 先实现 deterministic `analysis_story` builder；
5. 再扩展 plan 的 success criteria / stop conditions；
6. 最后实现 evidence-linked interpretation report。

完成 P1–P3 后，再重新评估是否进入 velocity / CellRank / spatial / multimodal。
