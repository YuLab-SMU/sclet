# sclet R-native AI 接入规划与实现 Spec

- **状态**：Phase 0–2 and the Phase 3 planning/controlled-execution foundation implemented; broader built-in action catalog remains future work
- **范围**：R 包内通过 `aisdk` 接入 AI，并将 AI 作为单细胞分析的审计、规划和解释层
- **第一阶段原则**：用户不离开 R；外部 agent 作为第二阶段客户端，不参与第一阶段核心设计
- **目标文件**：`R/ai-*.R`、`R/copilot.R`、`R/status.R`、`R/analysis-accessors.R`、`tests/testthat/`

---

## 1. 背景与问题定义

sclet 已经具备统一的 analysis-state 基础设施：

- `sclet_set_analysis_state()` 记录分析类型、方法、输入、参数、产物和摘要；
- `get_analysis_context()` 提供当前 active view 和部分分析记录；
- `Status()` 提供可用 assay、reduction、analysis、health 和 active 状态；
- 多条分析主线已经接入状态系统，包括 integration、annotation、trajectory、velocity、CellRank、SCENIC、gene-set scoring、perturbation、priority、rare cells 等；
- `R/copilot.R` 已有 `SummarizeContextForLLM()`、`sclet_copilot()` 和 `AuditAnalysisChain()` 的初步实现。

当前问题不是“没有 AI 入口”，而是 AI 接入仍然偏向一次性 prompt：

1. `SummarizeContextForLLM()` 将部分状态拼接为文本；
2. AI 主要只能生成解释，不能可靠地查询状态、检查前置条件或调用受控的 R 工具；
3. 账本记录了什么、结果由什么证据支撑、下一步允许做什么，没有统一的 AI contract；
4. `Status()` 已经能显示 `priority` / `rare_cells`，但 AI context 也必须同步纳入这些记录；
5. AI 输出缺少稳定的结构化结果，难以被后续 R 函数消费和写回 ledger；
6. AI 的分析建议、风险和结论边界没有形成可追溯记录。

本 Spec 的目标是将 AI 建设为 sclet 内部的一个 R-native 分析层：

> R 负责事实、统计计算、数据访问、状态更新和约束；AI 负责审计、规划、解释和自然语言交互。

---

## 2. 目标与非目标

### 2.1 目标

第一阶段完成以下能力：

- 在 R 内通过 `aisdk` 统一调用不同模型 provider；
- 为 AI 提供机器可读、版本化、可裁剪的 analysis ledger view；
- 将 AI 封装为职责清晰的 R 函数，而不是只提供一个泛化 chatbot；
- 为 AI 提供只读状态查询工具；
- 让 AI 输出结构化的 findings、evidence、warnings、recommendations 和 proposed actions；
- 支持分析计划生成、前置条件检查和 dry-run；
- 在用户明确确认后执行白名单中的 R 分析函数；
- 执行结果自动写回 analysis-state；
- 对 AI 结论保留证据、范围、模型和上下文版本；
- 保持无 AI 环境时 sclet 的所有确定性分析功能正常工作。

### 2.2 非目标

第一阶段不做：

- 不把外部 agent 作为主要入口；
- 不直接把整个 `SingleCellExperiment` 或表达矩阵序列化给模型；
- 不让模型自行计算关键统计量；
- 不让模型任意执行 R 代码或任意调用包函数；
- 不把模型生成的生物学解释自动当成事实；
- 不重写现有 analysis-state contract；
- 不在第一阶段实现多 agent 协作、长期向量记忆或复杂知识图谱；
- 不要求所有现有分析函数一次性补齐完整 provenance，先建立兼容层和最小字段。

---

## 3. 设计原则

### 3.1 Facts first

所有关键事实优先由 R 计算：

- 细胞数、基因数、assay/layer/reduction；
- 样本和分组设计；
- 统计检验、p 值、效应量和 QC 指标；
- 输入是否存在；
- 分析是否完成；
- 分析结果之间的差异和一致性。

AI 只解释 R 提供的事实，不负责替代统计计算。

### 3.2 Ledger is the source of truth

AI 不维护一份独立的事实状态。每次调用都从当前 SCE 和 analysis ledger 读取上下文；分析执行后，结果由 R 写回 ledger。

### 3.3 Structured context, natural-language answer

给模型的 context 使用结构化 list/JSON；给用户的结果可以包含自然语言，但底层必须保留结构化字段。

### 3.4 Read-only before write

先实现只读 AI 函数，再实现分析计划，最后实现需要确认的执行功能。默认不能由 AI 修改 SCE。

### 3.5 Evidence-bound claims

每个 AI finding 必须尽可能带：

- `evidence_refs`：支持它的分析、表格、图或统计摘要；
- `scope`：适用的样本、细胞群和数据层；
- `claim_level`：observed / measured / associated / consistent_with / causal；
- `uncertainty`：不确定性、缺失信息和替代解释。

### 3.6 Provider isolation

`aisdk` 的具体 model/provider API 只能出现在内部 adapter 中。上层函数不能依赖某一个 provider 的具体实现，以便更换 OpenAI、DeepSeek、Claude、本地模型或其他 aisdk provider。

---

## 4. 总体架构

```text
SingleCellExperiment
        |
        v
sclet analysis-state + deterministic summaries
        |
        v
GetAnalysisLedger() / sclet_ai_context()
        |
        +------------------------+
        |                        |
        v                        v
sclet_ai_tools()           sclet_ai_schema()
        |                        |
        +----------+-------------+
                   v
             sclet_ai_call()
                   |
                   v
                 aisdk
                   |
                   v
          sclet_ai_result (R object)
                   |
          +--------+---------+
          |                  |
          v                  v
      user/report       explicit record
                              |
                              v
                 analysis-state / AI record
```

外部 agent 的第二阶段接口：

```text
external agent -> stable R/API boundary -> GetAnalysisLedger / tools / plan
```

外部 agent 不直接读取 SCE 内部实现，也不依赖 `sclet_*` 私有对象结构。

---

## 5. Ledger AI View

### 5.1 用户接口

新增导出函数：

```r
GetAnalysisLedger <- function(
    object,
    detail = c("summary", "full"),
    target = NULL,
    include_artifacts = FALSE,
    include_data = FALSE
)
```

建议的内部函数：

```r
sclet_ai_context <- function(
    object,
    target = c(
        "status", "qc", "planning", "markers", "trajectory",
        "velocity", "perturbation", "report"
    ),
    detail = c("summary", "full")
)
```

`GetAnalysisLedger()` 是稳定的公共接口；`sclet_ai_context()` 负责根据任务裁剪 context。

### 5.2 最小返回结构

```r
list(
    schema_version = "1.0",
    dataset = list(
        n_cells = ncol(object),
        n_genes = nrow(object),
        assays = ..., 
        layers = ..., 
        reductions = ..., 
        modalities = ..., 
        columns = ...
    ),
    active_view = list(
        assay = ..., 
        layer = ..., 
        reduction = ..., 
        graph = ..., 
        ident = ...
    ),
    health = Status(object)$health,
    analyses = list(...),
    workflows = list(...),
    lineage = list(...),
    capabilities = list(...),
    blocked_actions = list(...),
    quality_checks = list(...),
    warnings = list(...)
)
```

### 5.3 analysis node 最小字段

现有记录字段保持兼容，同时建议统一到以下语义：

```r
list(
    id = "velocity_scvelo_1",
    type = "velocity",
    status = "completed",
    method = "velociraptor::scvelo",
    inputs = list(
        assay = "spliced",
        unspliced_assay = "unspliced",
        reduction = "PCA"
    ),
    parents = c("normalization_1", "pca_1"),
    params = list(...),
    artifacts = list(
        velocity = list(kind = "matrix", ref = "..."),
        embedding = list(kind = "reduction", ref = "UMAP")
    ),
    summary = list(...),
    validation = list(...),
    warnings = list(...),
    software = list(...),
    created_at = "..."
)
```

第一阶段不要求所有旧记录立刻包含全部字段。`GetAnalysisLedger()` 应提供 `NULL`、`unknown` 或兼容性标记，不伪造缺失信息。

### 5.4 priority / rare cells / workflow 纳入规则

AI view 必须纳入：

- `priority` state records；
- `rare_cells` state records；
- `state_priority` workflow record；
- `Status()$health$has_perturbation_priority`；
- `Status()$health$has_rare_cells`；
- `mainline_missing`；
- workflow artifacts 和 summary。

不能出现 `Status()` 已经能看到结果，但 AI context 看不到结果的分裂状态。

### 5.5 大对象策略

默认不向模型发送：

- 完整表达矩阵；
- 完整稀疏矩阵；
- 全量细胞级 metadata；
- 大型图片和长表格；
- 未裁剪的 `SingleCellExperiment`。

只发送：

- 统计摘要；
- 表格前 N 行和 schema；
- 分位数、分组汇总、QC 指标；
- artifact reference；
- 模型需要时通过工具按需获取的结果。

---

## 6. aisdk Adapter

### 6.1 内部接口

新增内部函数：

```r
sclet_ai_call <- function(
    task,
    context,
    tools = list(),
    schema = NULL,
    model = NULL,
    system_prompt = NULL,
    temperature = NULL,
    max_tokens = NULL,
    ...
)
```

职责：

1. 检查 `aisdk` 是否安装；
2. 解析默认 model 或用户传入 model；
3. 创建 agent/model；
4. 注入结构化 context；
5. 注册只读或执行工具；
6. 请求结构化输出；
7. 解析和验证结果；
8. 将 provider 错误转成 sclet 可识别的错误；
9. 记录 model、task、schema version 和耗时。

### 6.2 配置接口

建议支持：

```r
sclet_ai_config <- function(
    model = NULL,
    provider = NULL,
    api_key_env = NULL,
    default_temperature = 0,
    max_context_chars = 50000,
    max_tool_calls = 12,
    allow_execution = FALSE
)
```

配置不应把 API key 写入 SCE、analysis-state 或日志。

第一阶段可以复用 `OPENAI_MODEL` 和现有 aisdk 默认 model 逻辑；具体 provider 解析集中在 adapter 内。

### 6.3 错误处理

至少区分：

- `sclet_ai_missing_dependency`：未安装 `aisdk`；
- `sclet_ai_missing_model`：没有可用 model；
- `sclet_ai_provider_error`：provider/API 错误；
- `sclet_ai_invalid_output`：模型输出不符合 schema；
- `sclet_ai_context_error`：无法生成合法 context；
- `sclet_ai_tool_error`：工具参数或执行失败；
- `sclet_ai_execution_denied`：没有用户确认或不在白名单。

AI 不可用时，确定性 R 函数必须仍然可以运行。

---

### 6.4 aisdk 原生工具与结构化输出

实现依据：YuLab-SMU/aisdk 当前 README 的 tools/agents/structured output 章节，以及 `aisdk` 1.5.0 API（`tool()`, `create_agent()`, `z_*()`, `generate_object()`）。维护时以本地安装版本的函数签名和 capability 检测为准：<https://github.com/YuLab-SMU/aisdk>。


- `aisdk::tool()` 将 R 函数封装为模型可调用工具；
- `aisdk::create_agent(..., tools = ...)` 支持 agent tool loop；
- `aisdk::z_object()`、`z_array()`、`z_enum()` 等提供结构化 schema；
- `aisdk::generate_object(..., mode = "tool")` 支持通过 tool mode 获取结构化结果；
- `create_permission_hook()` 和 sandboxed execution 可作为后续执行权限层。

因此 sclet 的 adapter 应优先使用原生能力，而不是把工具仅作为 prompt 文本描述：

1. sclet 保持 provider-neutral 的内部 tool descriptor；
2. adapter 将 descriptor 转换为 `aisdk::tool()`；
3. `create_agent()` 接收转换后的 tools，并限制 `max_steps`；
4. 结构化结果优先使用 `z_*` schema 和 `generate_object(mode = "tool")`；
5. R 端仍保留最终结果校验，因为模型、provider 和 aisdk 版本可能存在能力差异；
6. 当运行环境中的 aisdk 缺少某项能力时，adapter 可以退化到文本响应，但必须在 result metadata 中标记 `native_tools = FALSE` 或 `structured_output = FALSE`。

实现不能假设仅有 OpenAI；模型应使用 aisdk 的 provider:model 标识和默认 model 配置。

### 6.5 双阶段调用协议

当一个任务需要查询 ledger/tool 后再形成可持久化的结构化结论时，使用双阶段协议：

```text
Phase A: create_agent + native tools + bounded max_steps
        ↓
        tool-loop response
        ↓
Phase B: generate_object(mode = "tool") + z_* schema
        ↓
        R 端 sclet_ai_result 校验
```

`sclet_ai_call()` 的 `structured_output = TRUE` 启用该协议。第一阶段负责查询和推理，第二阶段只负责把第一阶段的回答压缩成结构化的 `answer/findings/evidence/warnings/recommendations/proposed_actions`。第二阶段失败时，默认保留第一阶段文本结果，并在 `metadata$structured_output = FALSE` 和 warnings 中明确记录 fallback；需要严格结构化时可设置 `fallback_on_structure_error = FALSE`。

高层只读函数默认启用结构化输出，但仍允许调用者显式关闭：

```r
AIStatus(sce, structured_output = FALSE)
```

这一步不执行新的单细胞分析，也不改变 SCE；它只是让 AI 结果可以可靠地被 `RecordAIResult()` 写回 ledger。

---

## 7. 结构化 AI 结果

### 7.1 结果类

新增内部构造函数：

```r
new_sclet_ai_result <- function(
    task,
    answer = NULL,
    findings = list(),
    evidence = list(),
    warnings = list(),
    recommendations = list(),
    proposed_actions = list(),
    context = NULL,
    metadata = list()
)
```

返回对象：

```r
class(result)
# c("sclet_ai_result", "list")
```

### 7.2 统一字段

```r
list(
    task = "qc_review",
    answer = "...",
    findings = list(
        list(
            id = "finding_1",
            severity = "warning",
            statement = "...",
            evidence_refs = c("qc_1", "status_1"),
            scope = list(...),
            claim_level = "observed",
            uncertainty = "..."
        )
    ),
    evidence = list(...),
    warnings = list(...),
    recommendations = list(...),
    proposed_actions = list(...),
    context = list(
        schema_version = "1.0",
        fingerprint = "..."
    ),
    metadata = list(
        model = "...",
        provider = "...",
        created_at = "...",
        duration_sec = ...
    )
)
```

### 7.3 输出验证

AI 输出必须先通过 R 端校验：

- 必需字段存在；
- `severity`、`claim_level` 使用允许值；
- `evidence_refs` 指向当前 context 中存在的对象，或明确标记为 unresolved；
- proposed action 的 action name 在 registry 中；
- 不能把 `causal` 作为默认 claim level；
- 缺失证据时必须降低结论等级或输出 warning。

---

## 8. 专用 AI 函数

第一阶段建议实现以下函数，全部默认为只读：

### 8.1 状态与 QC

```r
AIStatus(object, model = NULL, ...)
AIReviewQC(object, model = NULL, ...)
```

`AIReviewQC()` 输入由 R 预先生成 QC summary，AI 负责归纳问题、排序风险和提出建议。

### 8.2 分析解释

```r
AIExplainAnalysis(object, type = NULL, id = NULL, model = NULL, ...)
AICompareAnalyses(object, x, y, model = NULL, ...)
AIAuditAnalysisChain(object, target, model = NULL, ...)
```

现有 `AuditAnalysisChain()` 保留兼容，但内部应逐步迁移到统一 `sclet_ai_call()` 和结构化结果。

### 8.3 下一步和规划

```r
AIRecommendNextStep(object, model = NULL, ...)
AIPlanAnalysis(object, question, model = NULL, ...)
ValidateAIPlan(object, plan, ...)
```

计划至少包含：

```r
list(
    goal = "...",
    assumptions = list(...),
    steps = list(...),
    prerequisites = list(...),
    risks = list(...),
    expected_outputs = list(...),
    requires_confirmation = TRUE
)
```

### 8.4 领域解释

可在基础设施稳定后添加：

```r
AIInterpretMarkers(object, markers, group_by = NULL, ...)
AIInterpretTrajectory(object, id = NULL, ...)
AIInterpretPerturbation(object, id = NULL, ...)
AIInterpretCommunication(object, id = NULL, ...)
AIInterpretGeneSets(object, id = NULL, ...)
```

这些函数不应重新计算核心结果，而是调用 R 端已经生成的统计摘要、结果表和 evidence refs。

---

## 9. Tool Registry

### 9.1 第一阶段只读工具

建议建立内部 registry：

```r
sclet_ai_tool_registry <- function(object, mode = c("read", "execute"))
```

第一批只读工具：

```text
get_status
get_ledger
get_capabilities
get_analysis_record
get_analysis_lineage
get_analysis_dependents
get_quality_checks
get_evidence_summary
compare_analysis_states
fetch_artifact_summary
```

每个工具必须定义：

```r
list(
    name = "get_analysis_record",
    description = "...",
    input_schema = list(...),
    output_schema = list(...),
    read_only = TRUE,
    max_output_size = ...,
    handler = function(...) ...
)
```

### 9.2 执行工具

第二阶段才允许注册：

```text
RunVelocity
RunTrajectoryWorkflow
RunStatePriorityWorkflow
RunGeneSetScoring
RunMilo
RunMarkers
```

执行工具必须满足：

- 白名单；
- 明确 input schema；
- 前置条件检查；
- dry-run；
- 用户确认；
- 参数可追溯；
- 出错时不写入伪造的 completed record；
- 成功后自动注册 analysis-state；
- 最好通过 plan hash 防止执行过期计划。

禁止：

- `eval(parse(...))`；
- AI 传入任意 R 表达式；
- AI 任意指定文件路径覆盖用户文件；
- 未确认时直接运行昂贵或有副作用的分析。

---

## 10. 计划与执行协议

### 10.1 规划阶段

```r
plan <- AIPlanAnalysis(
    sce,
    question = "比较不同 condition 下的细胞状态变化"
)
```

规划输出必须能被 R 验证，而不是只有自然语言：

```r
validation <- ValidateAIPlan(sce, plan)
```

返回：

```r
list(
    valid = FALSE,
    steps = list(...),
    missing = list("sample_id", "condition"),
    warnings = list("当前重复数未确认"),
    estimated_cost = list(...),
    requires_confirmation = TRUE
)
```

### 10.2 执行阶段

建议接口：

```r
sce2 <- ExecuteAIPlan(
    sce,
    plan = plan,
    confirmation = TRUE,
    dry_run = FALSE
)
```

第一版可以要求用户显式设置：

```r
ExecuteAIPlan(sce, plan, confirm = TRUE)
```

执行过程：

1. 冻结当前 context fingerprint；
2. 再次验证 prerequisites；
3. 检查 plan 中所有 action 是否在 registry；
4. 按步骤调用 R 函数；
5. 每一步成功后写入 analysis-state；
6. 失败则记录失败节点和错误，不伪装成完成；
7. 返回更新后的 SCE 和 execution report。

---

## 11. AI 结果写回账本

AI 分析默认不改变对象。需要明确写回时使用：

```r
RecordAIResult <- function(object, result, id = NULL, active = FALSE)
```

记录类型建议：

```text
ai_review
ai_plan
ai_interpretation
ai_audit
ai_execution
```

记录内容至少包含：

```r
list(
    task = result$task,
    model = result$metadata$model,
    context_schema_version = result$context$schema_version,
    context_fingerprint = result$context$fingerprint,
    evidence_refs = ...,
    findings = ...,
    recommendations = ...,
    proposed_actions = ...,
    created_at = ...
)
```

不能把 API key、完整 prompt 中的敏感信息或无法复现的临时对象写入 ledger。

---

## 12. 科研主张约束

AI 输出必须区分以下层级：

```text
observed       观察到
measured       测量/统计得到
estimated      估计
associated     相关
consistent_with 与某机制一致
suggestive     提示性证据
causal         因果性，仅在有明确设计和证据时允许
```

默认规则：

- 普通单细胞表达差异默认不升级为 `causal`；
- trajectory 不自动等同于真实时间因果过程；
- velocity 不自动证明细胞命运改变；
- integration 后的 embedding 不自动等同于 corrected expression；
- marker/pathway 富集结果必须带数据范围和背景集；
- 没有样本级重复信息时，AI 必须明确降低推断强度；
- 没有 evidence ref 的生物学结论应标记为 hypothesis，而非 finding。

---

## 13. 测试策略

### 13.1 无模型单元测试

不请求真实模型，测试：

- `GetAnalysisLedger()` 的 schema；
- priority / rare cells / state_priority 是否进入 AI view；
- 空对象和部分状态对象的兼容性；
- context detail 和 target 裁剪；
- artifact 大小限制；
- lineage 和 evidence refs；
- plan validation；
- tool registry schema；
- AI result 构造和字段校验；
- claim level 和 evidence 约束；
- AI record 写回 analysis-state。

### 13.2 aisdk mock 测试

使用 mock adapter，不调用网络：

- 模型返回合法结构化结果；
- 模型返回缺字段结果；
- 模型返回非法 action；
- provider 错误；
- tool 调用错误；
- 超出最大 context 或最大 tool call 数；
- 没有 aisdk 时错误类型正确。

### 13.3 执行测试

使用轻量或 mock R analysis function：

- plan 中存在缺失 prerequisite 时不执行；
- 未确认时不执行；
- 非白名单 action 被拒绝；
- 中途失败时不产生 completed record；
- 成功后 analysis-state 正确更新；
- 重复执行不会破坏已有记录；
- context fingerprint 变化后旧计划被拒绝或要求重新验证。

### 13.4 推荐测试文件

```text
tests/testthat/test-ai-context.R
tests/testthat/test-ai-result.R
tests/testthat/test-ai-tools.R
tests/testthat/test-ai-planning.R
tests/testthat/test-ai-execution.R
tests/testthat/test-ai-copilot.R
```

现有 `test-state-refactor.R` 继续覆盖 `Status()`；不把全部 AI 行为塞进该文件。

---

## 14. 分阶段实现计划

### Phase 0：契约和兼容层

目标：先让 AI 能完整、稳定地看到 ledger。

任务：

1. 新增 `GetAnalysisLedger()`；
2. 新增 `sclet_ai_context()`；
3. 纳入所有已有 analysis records；
4. 纳入 `state_priority`、`priority`、`rare_cells` 和 `Status()$health`；
5. 为旧记录提供兼容 normalization；
6. 将 `SummarizeContextForLLM()` 改为基于 AI context 生成文本；
7. 加 schema 和 fingerprint。

验收：

- 同一个 SCE 在 summary/full view 下结果稳定；
- AI view 不丢失 `Status()` 已显示的状态；
- 不发送完整表达矩阵；
- 无模型时 context 函数仍然可用。

### Phase 1：aisdk adapter 和只读 AI

任务：

1. 新增 `sclet_ai_call()`；
2. 新增 `sclet_ai_result`；
3. 新增 schema 校验；
4. 实现 `AIStatus()`；
5. 实现 `AIReviewQC()`；
6. 实现 `AIExplainAnalysis()`；
7. 将 `sclet_copilot()` 迁移到统一 adapter。

验收：

- 上层 API 不依赖具体 provider；
- mock 模型可测试；
- 输出同时包含自然语言和结构化 findings/evidence/warnings。

### Phase 2：工具查询和审计

任务：

1. 建立 read-only tool registry；
2. 实现 lineage、evidence、capabilities 查询；
3. 实现 `AIAuditAnalysisChain()`；
4. 实现 `AIRecommendNextStep()`；
5. 实现 `RecordAIResult()`。

验收：

- AI 可以按需查询，而非依赖单个超长 prompt；
- 每个建议都能关联到当前 ledger；
- AI audit 可回写为可追溯记录。

### Phase 2.5：双阶段结构化总结（已实现）

在只读 tool registry 基础上，`sclet_ai_call(structured_output = TRUE)` 先运行 aisdk native tool loop，再使用 `generate_object(mode = "tool")` 和 `z_*` schema 将最终回答归一化为 `sclet_ai_result`。高层只读 AI 函数默认启用该模式；结构化阶段失败时保留第一阶段文本并记录 fallback warning，除非调用者设置 `fallback_on_structure_error = FALSE`。

验收：

- native tool descriptors 能转换为 aisdk `Tool` 对象；
- mock agent + mock structured endpoint 能覆盖双阶段协议；
- 结构化结果由 R 端再次校验；
- fallback 不伪造结构化成功；
- `make check` 通过且无 errors/warnings/notes。

---

## Phase 3：分析规划和受控执行（基础版已实现）

任务：

1. 实现 `AIPlanAnalysis()`；
2. 实现 `ValidateAIPlan()`；
3. 建立 execution action registry；
4. 实现 dry-run 和确认机制；
5. 实现 `ExecuteAIPlan()`；
6. 自动登记执行结果和失败记录。

验收：

- AI 不能调用未注册函数；
- 缺失前置条件时不会执行；
- 用户确认前不会执行有副作用操作；
- 每个结果都回到账本。

Phase 3 基础实现采用显式 allowlist，而不是从全局环境按字符串查找函数：

- `AIAction()` 描述一个可执行 handler、参数、前置条件和返回类型；
- `AIExecutionRegistry()` 只接受显式注册的 action，默认是空 registry；
- `AIDefaultExecutionRegistry(object)` 提供 `inspect_status`、`inspect_ledger`、`check_qc` 三个只读 action，它们不修改 SCE，也不要求 confirmation；
- `AIAction()` 还记录 input/output schema、是否修改对象、允许写入的 state 类型、估计成本和幂等性；
- `AIPlanAnalysis()` 只生成 `sclet_ai_plan`，不会执行任何 action；
- `ValidateAIPlan()` 检查 action 是否注册、参数 schema、步骤依赖、当前 ledger fingerprint 和前置条件；
- `ExecuteAIPlan()` 默认 `dry_run = TRUE`，包含修改对象或要求确认的 action 时，非 dry-run 必须提供 validation 返回的 confirmation token；
- 成功或失败执行都可以通过 `sclet_set_analysis()` 写入 `ai_execution_*` 记录，并写入 command log；
- execution registry 不作为 AI tool 暴露，AI 不能自行触发有副作用的函数。

Phase 3.1 暂不提供默认有副作用 action。具体分析函数必须由用户或上层 workflow 通过 `AIAction()` 显式包装并注册，后续可以在稳定的 registry contract 上增加经过审查的 sclet-native actions。

### Phase 4：外部 agent 接入

- 以 `GetAnalysisLedger()` 作为外部 context boundary；
- 以 tool registry 作为外部工具描述；
- 以 plan validation / execution protocol 作为安全边界；
- 外部 agent 不读取 SCE 内部 slots，也不绕过 R 端验证。

---

## 15. 建议的文件组织

```text
R/
  ai-context.R       # GetAnalysisLedger, sclet_ai_context
  ai-adapter.R       # sclet_ai_call, aisdk/provider adapter
  ai-result.R        # sclet_ai_result and validation
  ai-tools.R         # read-only and execution tool registry
  ai-planning.R      # AIPlanAnalysis, ValidateAIPlan
  ai-execution.R     # ExecuteAIPlan, RecordAIResult
  ai-functions.R     # AIStatus, AIReviewQC, AIExplainAnalysis, ...
  copilot.R          # backward-compatible public facade

tests/testthat/
  test-ai-context.R
  test-ai-result.R
  test-ai-tools.R
  test-ai-planning.R
  test-ai-execution.R
  test-ai-copilot.R
```

不建议第一阶段把所有内容继续堆进 `copilot.R`。保留该文件作为兼容入口，内部实现迁移到拆分后的模块。

---

## 16. 第一版完成标准

以下条件全部满足，才认为 R-native AI 第一版完成：

- [ ] `GetAnalysisLedger()` 能输出版本化结构化 context；
- [ ] `priority`、`rare_cells`、`state_priority` 和 health 信息对 AI 可见；
- [ ] `SummarizeContextForLLM()` 基于统一 context 生成；
- [ ] `sclet_ai_call()` 隔离 aisdk provider 细节；
- [ ] `sclet_ai_result` 有稳定字段和本地校验；
- [ ] 至少实现 `AIStatus()`、`AIReviewQC()`、`AIExplainAnalysis()`、`AIRecommendNextStep()`；
- [ ] 有 mock 测试，不依赖网络和真实模型；
- [ ] AI 输出带 evidence、scope、claim level 或明确的缺失标记；
- [ ] 默认只读，不会未经确认修改 SCE；
- [ ] 重要 AI 结果可以显式写回 analysis-state；
- [ ] 没有 `aisdk` 时，普通 sclet 分析和状态接口不受影响；
- [ ] 外部 agent 可以在不理解 SCE 内部结构的情况下，消费稳定的 ledger/tool contract。

---

## 17. 设计结论

sclet 的 AI 不应只是一个“把状态拼到 prompt 后问模型”的 Copilot，而应成为一个 R-native 的分析协作者：

1. 通过 ledger 了解当前事实和血缘；
2. 通过 tools 查询证据和前置条件；
3. 通过专用函数完成审计、规划和解释；
4. 通过 plan/validate/execute 逐步进入受控分析执行；
5. 通过 analysis-state 持久化 AI 工作痕迹；
6. 未来再把同一套 contract 暴露给外部 agent。

第一步不是增加更多模型，而是先把 `GetAnalysisLedger()`、`sclet_ai_context()`、`sclet_ai_call()` 和结构化 `sclet_ai_result` 建起来。它们将成为后续所有 R-native AI 函数和外部 agent 接入的共同基础。
