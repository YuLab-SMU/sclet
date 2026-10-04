# sclet AI 高级个性化分析开发 Spec

- **状态**：Phase C integration 真实多路线执行、comparison trade-off 检测、evidence dependency/scope 收口、privacy 默认严格门控 已实现；Phase D 第一批 marker/DE/annotation 证据链（`check_annotation_readiness()`、`run_de_test`、`run_annotation`）已实现；Phase D 第二批 rare-cell/doublet 证据链（`check_rare_cell_readiness()`、`summarize_small_cluster_evidence()`、`rare_cell` action group）与 Phase D 第三批 trajectory 第一批（`check_trajectory_readiness()`、只读 `summarize_trajectory_cluster_order()`、`trajectory` action group 与 root 必填的 `run_trajectory`）已实现；velocity / CellRank / fate / spatial / multimodal action catalog、trajectory 多 root 比较、完整 causal ceiling 仍在规划中
- **版本**：0.3
- **范围**：围绕 `SingleCellExperiment`、analysis ledger 和 `aisdk`，建设面向高级单细胞问题的 AI 主导分析层
- **核心判断**：标准化基础分析优先使用确定性的 R workflow；AI 的主要价值放在数据诊断、个性化路线选择、高级分析编排和证据整合
- **当前基础**：`GetAnalysisLedger()`、结构化 `sclet_ai_result`、`AIPlanAnalysis()`、`ValidateAIPlan()`、`ExecuteAIPlan()`、`RunAIPlan()`、`RunAIAnalysis()`、`AskAI()` 和 native action registry 已实现
- **不在本 Spec 内**：外部 agent、多 agent 协作、任意 R 代码执行、把完整表达矩阵直接发送给模型
- **重要状态说明**：本文是高级 AI 分析架构提案与分阶段路线图。`GetAIProfile()`、第一批 deterministic diagnostics、`AIInvestigate()`、`RecordAIEvidence()`、`ValidateAIEvidenceRefs()`、payload sanitizer、只读 `CompareAIAnalyses()`、`RunIntegrationRoutes()`（真实多路线执行 + baseline/metrics 登记）、`AIDefaultExecutionRegistry` integration/annotation group、`ConfirmAIDesignSemantics()`（骨架指纹 + 逐列值签名双重防伪的 design 确认闸门）、`check_annotation_readiness()`、`run_de_test`/`run_annotation`（reference 显式必填、evidence claim_level 限定为 associated/consistent_with、不覆盖 Idents、不含逐细胞原始值）、`sclet_ai_evidence_independence()`、默认 `enforce_privacy=TRUE` 的 `sclet_ai_call` UX、`dependency_group` hash、scope 越界与 cross-kind 互引拒绝 均已实现；velocity / CellRank / fate / spatial / multimodal action catalog、trajectory 多 root 比较、完整因果 claim ceiling、多 agent 协作仍未实现。
- **trajectory 已实现范围（仅第一批）**：`check_trajectory_readiness()`（只回答“cluster + 可用 embedding 是否齐备”，**不输出任何起点/root 建议**）、只读的 `summarize_trajectory_cluster_order()`（只描述 cluster 在嵌入维度上的分布，`root_suggested` 恒为 `FALSE`）、`AIDefaultExecutionRegistry` 的 `trajectory` group（`run_trajectory`，包装 `RunSlingshot()`，`group` 与 `start_cluster` 均必填，缺失/NULL/空/不存在的 root 分别被 `start_cluster_missing`/`start_cluster_unknown` 拒绝）、以及 `claim_level = "consistent_with"` 的聚合 evidence（只含 lineage 数量、pseudotime 分位数与匿名化 `cluster_N` 起点编码，不含逐细胞向量，`pseudotime_is_absolute_time = FALSE`）。**尚未覆盖**：`compare_trajectory_roots` 多 root 比较、`check_velocity_readiness`、`run_velocity`、`run_fate_analysis`、spatial 和 multimodal 的任何执行 action（需要 spliced/unspliced 或 Python 后端，留给后续批次）。
- **rare-cell / doublet 证据链已实现范围**：`check_rare_cell_readiness()`（cluster 分配、PCA reduction、doublet 证据缺失的显式说明）、只聚合不下结论的 `summarize_small_cluster_evidence()`（QC / doublet / marker / sample replication 四类独立信号 + 计数）、`AIDefaultExecutionRegistry` 的 `rare_cell` group（`run_doublet_detection`、`run_rare_cell_detection`）、按独立信号数量分级的 evidence `claim_level`（0 个信号不登记 evidence、1 个信号降级为 `associated` 并标 `low_confidence = TRUE`、≥2 个信号才 `consistent_with`），以及“只标注不删除”的红线。**尚未覆盖**：ambient RNA（decontX）信号本轮刻意不纳入独立证据，因为 `RunDecontX()` 的 state 登记尚未加固；`compare_rare_cell_evidence` 尚未实现（多路线 rare-cell 比较）；删除/合并稀有群体的独立 action 明确不在本轮范围内。文中“建议”“拟支持”“应”表示未来 contract；“当前”只指本文件列出的已实现能力。

## 0. 本轮审核结论与采纳范围

本轮审核反馈大部分采纳，主要修订如下：

- **采纳**：把高级能力明确标为未来能力，避免将 proposed API 写成已实现 API；
- **采纳**：增加 outbound payload 和字段级隐私 contract；
- **采纳**：增加 evidence ref 验证、claim ceiling、scope 和证据独立性；
- **采纳**：增加 `clarification_required`，禁止 AI 猜测 batch/condition/subject 语义；
- **采纳**：明确执行状态、部分失败和完整事务回滚的差异；
- **采纳**：增加 integration、annotation、rare-cell 的可验收指标 contract；
- **采纳**：补充 privacy、stale plan、过期确认、evidence 和 partial failure 测试；
- **采纳**：将第一批开发任务收窄为 profile、诊断、integration readiness、调查和路线比较；
- **暂不承诺**：完整事务回滚、causal claim、全量高级 action catalog 和外部 agent；这些保留为后续研究或工程阶段。

---

## 1. 背景与产品定位

### 1.1 基础分析不是主要 AI 场景

以下流程已经有成熟的确定性实现：

```text
QC → normalization → HVG → scaling → PCA → neighbors → clustering → UMAP
```

它的特点是步骤相对固定、参数有成熟默认值、错误主要属于工程性前置条件错误。AI 可以帮助新手理解这些步骤，但不应为了体现 AI 而重新规划一条已经稳定的标准路径。

已提供基础 workflow facade（`RunBasicWorkflow()`，纯确定性 R 组合，不涉及 AI）：

```r
sce <- RunBasicWorkflow(sce)  # 已实现，可选 steps = / n_pcs / cluster_resolution
```

该接口按标准顺序串联 `NormalizeData()` → `FindVariableFeatures()` → `ScaleData()` → `RunPCA()`
→ `FindNeighbors()` → `FindClusters()` → `RunUMAP()`，每一步都调用现有函数、由其自行登记 state，
facade 本身不重复登记、不重新实现任何算法、不进入 `AIDefaultExecutionRegistry()`。
专家用户仍应直接调用底层函数。

AI 在此阶段主要用于解释状态、解释图形、提示明显问题，而不是作为主要决策者。

### 1.2 AI 的核心价值

AI 的主要价值应放在以下问题：

- 当前数据最值得优先解决的问题是什么？
- 当前 cluster 是真实生物群体、批次效应、doublet 还是低质量细胞？
- 是否应该做 integration，应该选择哪一种策略？
- 某个小群体是否值得保留和进一步验证？
- 聚类、marker、样本分布和实验设计是否相互支持？
- 当前结果支持 annotation、trajectory、DE 或 communication 吗？
- 多条分析路线的证据是否一致？
- 下一步分析的收益、成本和风险分别是什么？

因此产品定位应从：

> AI 帮新手运行标准单细胞函数

转向：

> AI 读取数据状态和分析证据，为每个数据集选择并编排个性化的高级分析路线；R/sclet 负责计算、验证和记录，人负责科学目标与关键边界。

### 1.3 三层用户模式

| 模式 | 入口 | 主要用户 | AI 角色 |
|---|---|---|---|
| Basic | `RunBasicWorkflow()` | 新手和标准项目 | 解释为主 |
| Guided | `RunAIAnalysis()`、`AskAI()` | 新手和一般用户 | 诊断、建议、带确认执行 |
| Advanced | `AIInvestigate()`、高级 plan/action API | 研究人员 | 个性化路线设计、多方案比较、证据整合 |

---

## 2. 设计目标与非目标

### 2.1 目标

1. AI 可以识别数据集的实验结构、模态、样本分组和当前分析状态。
2. AI 可以根据数据特征而不是固定模板提出高级分析路线。
3. 每个高级建议都能关联到 R 计算出的事实、统计摘要或已有分析记录。
4. 高级分析 action 使用可验证的输入 schema、输出 contract 和 prerequisite。
5. 多个分析分支可以并行或顺序执行，并且能够比较结果。
6. 分析结果、AI 判断、用户决策和替代方案全部进入可追溯 ledger。
7. 新手不需要直接面对 registry、fingerprint 和 confirmation token。
8. 专家可以访问完整 plan、validation、execution 和 evidence graph。
9. 没有 API key 时，确定性分析和离线测试仍然完全可用。
10. 默认模型由 aisdk 配置；上层函数不要求用户重复传 `model`。

### 2.2 非目标

- 不让 AI 代替统计方法或直接计算 p 值、效应量、显著性和距离矩阵。
- 不允许 AI 生成并执行任意 R/Python 代码。
- 不把完整 counts、表达矩阵或受保护的样本信息默认放进模型上下文。
- 不将 AI 的 cell-type 命名自动写成事实。
- 不在本阶段实现外部 agent 或多 agent 协作。
- 不试图一次性覆盖所有 Bioconductor 和第三方单细胞包。
- 不把所有基础函数重新包一层而缺少高级诊断价值。

---

## 3. 目标用户体验

### 3.1 基础 workflow

标准分析未来由确定性 workflow facade 完成：

```r
sce <- RunBasicWorkflow(sce)  # 拟议接口，当前版本尚未实现
```

当前版本仍使用已有的确定性函数完成基础流程；本 Spec 不把该 facade 当作已实现 API。

AI 只在需要时解释：

```r
AskAI(sce, "请解释当前 PCA、cluster 和 UMAP，并指出明显风险。")
```

### 3.2 高级个性化问题

用户提供科学问题，而不是函数名称：

```r
result <- RunAIAnalysis(
    sce,
    goal = paste(
        "请判断这批数据是否存在批次效应。",
        "如果存在，请比较适合的 integration 路线，",
        "并告诉我每条路线需要什么证据来验证。"
    )
)
```

用户看到的交互应该是：

```text
AI 发现：PCA 前 3 个成分与 sample_id 的关联高于与当前 cluster 的关联。

可能原因：
1. 批次效应；
2. 样本组成差异；
3. 真实的条件特异性生物效应。

建议路线：
A. Harmony：成本低，主要校正 embedding；
B. fastMNN：适合校正低维结构；
C. scVI：适合更复杂的非线性批次结构，但计算成本较高。

需要你确认：
- sample_id 是否为技术批次？
- condition 是否为主要生物变量？
- 是否允许校正 condition？
```

AI 不应直接跳到“执行 Harmony”。它需要先提出问题、说明替代解释，并请求实验设计信息。

### 3.3 用户最小输入

高级分析至少需要：

```r
list(
    object = sce,
    goal = "自然语言科学问题",
    design = list(
        sample = "sample_id",
        batch = "batch_id",
        condition = "condition",
        subject = "patient_id"
    ),
    constraints = list(
        preserve_condition = TRUE,
        max_runtime = "medium",
        no_python = FALSE
    )
)
```

如果 `design` 缺失，AI 必须先询问，而不能自行猜测 `batch`、`condition` 或 `subject`。

### 3.4 当前实现与未来 API 的边界

当前已经可以使用：

```r
RunAIAnalysis(sce, goal = "...", confirm = "ask")
AskAI(sce, "...")
AIPlanAnalysis(sce, goal = "...")
ValidateAIPlan(plan, object = sce, registry = registry)
ExecuteAIPlan(...)
```

当前已经实现：

```r
GetAIProfile(sce)
AIInvestigate(sce, question = "...")
summarize_qc_by_group(sce, group = "sample_id")
check_integration_readiness(sce, design = list(...))
```

当前仍未实现、仅在本 Spec 中提出的 API 包括：

```r
CompareAIAnalyses()
RunBasicWorkflow()
```

当前版本的 `AIInvestigate()` 是只读骨架：它可以组合 bounded profile、ledger 和 readiness diagnostics，但不会执行高级 action，也不会自动生成完整的多路线比较。

当前包中已有部分高级底层算法，例如 integration、marker/DE、rare-cell、trajectory、velocity、spatial 和 multimodal 相关函数。后续重点不是重复实现这些算法，而是为它们补充 AI-native adapter、设计语义校验、确定性诊断、结果比较、evidence refs 和 provenance contract。

---

## 4. 总体架构

```text
SingleCellExperiment
        |
        v
analysis-state + deterministic summaries
        |
        +-----------------------------+
        |                             |
        v                             v
Dataset Profiler                GetAnalysisLedger
        |                             |
        +-------------+---------------+
                      v
              AI Investigation Context
                      |
          +-----------+------------+
          |                        |
          v                        v
   Read-only tools          Advanced action registry
          |                        |
          +-----------+------------+
                      v
                aisdk adapter
                      |
                      v
             structured AI result
                      |
       +--------------+---------------+
       |                              |
       v                              v
  user explanation             validated plan
                                      |
                              dry-run + confirmation
                                      |
                                      v
                              R/sclet execution
                                      |
                                      v
                              ledger + evidence graph
```

### 4.1 分层职责

#### A. Dataset profiler

由 R 确定性计算：

- 细胞、基因、样本和模态数量；
- assays、layers、reducedDims、graphs；
- `colData` / `rowData` 字段类型和缺失率；
- sample、batch、condition、subject 候选字段；
- counts 稀疏度和矩阵规模；
- QC 指标可用性；
- 已有分析和未满足的主线；
- 可能的分析成本。

#### B. AI investigation layer

AI 负责：

- 选择要检查的 deterministic tool；
- 组合事实形成问题诊断；
- 识别多种可能解释；
- 生成候选路线；
- 明确不确定性和缺少的信息；
- 请求用户补充实验设计。

#### C. Advanced action layer

action 负责真实计算，例如：

- batch diagnostics；
- integration；
- marker / DE；
- annotation；
- rare-cell detection；
- trajectory；
- velocity；
- pathway / gene program；
- cell-cell communication；
- spatial niche；
- multimodal analysis。

#### D. Evidence and ledger layer

所有 action 输出都必须生成：

- state record；
- command record；
- input references；
- output references；
- summary statistics；
- warnings；
- evidence references；
- software and parameter metadata。

## 5. 数据隐私与 outbound payload contract

高级 AI 分析必须先定义模型请求允许携带什么。仅声明“不发送完整矩阵”不够；需要字段级策略。

### 5.1 默认允许发送

默认只允许发送：

- 细胞、基因、样本数量的聚合摘要；
- assay/layer/reduction/graph 名称；
- 经过聚合的 QC、PCA、cluster 和样本组成统计；
- 软件版本、参数和状态记录；
- 已脱敏的 group labels（例如 `group_1`）；
- 用户明确提供的自然语言目标。

### 5.2 默认禁止发送

默认禁止：

- 完整 counts 或表达矩阵；
- cell barcode；
- patient/subject/sample 原始 ID；
- 自由文本形式的临床信息；
- 小于隐私阈值的可识别群体计数；
- 未经用户确认的 `colData` 原始值；
- analysis record 中可能包含的 prompt、路径、token 或密钥。

### 5.3 payload policy

新增内部策略函数：

```r
sclet_ai_payload_policy <- function(
    object,
    context,
    privacy = c("strict", "standard", "local")
)
```

策略必须返回：

```r
list(
    allowed_fields = c(...),
    redacted_fields = c(...),
    aggregation_threshold = 10L,
    estimated_tokens = ..., 
    requires_user_consent = FALSE,
    warnings = c(...),
    payload_fingerprint = "..."
)
```

### 5.4 安全验收

- payload 发送前必须经过 policy 检查；
- `AskAI()` 和高级调查函数不能绕过 policy；
- metadata 中的自然语言内容必须标记为 data，不得当作系统指令；
- payload 中不能出现 `DEEPSEEK_API_KEY`、provider token 或本地私密路径；
- strict 模式下发现未分类字段时拒绝发送，而不是静默放行；
- mock 测试必须能检查最终 outbound context，而不是只检查上层参数。

---

## 6. AI Context 扩展

当前 `GetAnalysisLedger()` 已经提供基础的状态视图。高级个性化分析需要增加一个面向诊断的 context 层，但仍然不能包含完整表达矩阵。

建议新增：

```r
GetAIProfile <- function(
    object,
    design = NULL,
    include_diagnostics = TRUE,
    include_cost = TRUE
)
```

### 6.1 profile 最小结构

```r
list(
    schema_version = "2.0",
    dataset = list(
        n_cells = ..., 
        n_features = ..., 
        assays = ..., 
        layers = ..., 
        reductions = ..., 
        modalities = ...
    ),
    design = list(
        sample_candidates = ..., 
        batch_candidates = ..., 
        condition_candidates = ..., 
        subject_candidates = ..., 
        missing_fields = ...
    ),
    qc = list(
        available = ..., 
        summaries = ..., 
        thresholds_used = ..., 
        filtering_history = ...
    ),
    structure = list(
        sample_sizes = ..., 
        group_sizes = ..., 
        sparsity = ..., 
        feature_detection = ...
    ),
    analysis_state = GetAnalysisLedger(object),
    diagnostics = list(...),
    capabilities = list(...),
    cost_estimates = list(...),
    privacy = list(...)
)
```

### 6.2 deterministic diagnostics

第一版只加入低成本、有明确解释的诊断：

- sample / batch / condition 的细胞数量；
- QC 指标在样本间的分布摘要；
- top PCA components 与 metadata 的关联摘要；
- cluster 与 sample 的交叉分布；
- cluster size、singleton 和极小群体比例；
- marker overlap / contradiction 摘要；
- embedding 前后或校正前后的结构差异；
- integration 前后的 batch mixing 和 biological separation 摘要；
- trajectory 输入是否完整；
- RNA velocity 所需 spliced / unspliced 是否存在。

模型只解释这些摘要，不自行从矩阵猜结论。

---

## 7. 高级 action catalog 路线

高级 action 不一次性全部实现，按科研价值和可验证性分阶段建设。

### P0：数据诊断与路线选择

目标：让 AI 能判断“当前数据最需要检查什么”。

建议 action：

```text
profile_dataset
inspect_design
summarize_qc_by_group
summarize_pca_metadata_association
summarize_cluster_sample_composition
summarize_small_clusters
check_integration_readiness
check_annotation_readiness
check_trajectory_readiness
check_velocity_readiness
```

这些 action 默认只读，不修改对象。

### P1：批次、整合与结果稳定性

建议 action：

```text
run_batch_diagnostic
run_integration
compare_integrations
compare_pre_post_correction
assess_cluster_stability
```

关键约束：

- 必须区分技术 batch 和生物 condition；
- 默认禁止校正用户标记为必须保留的 biological variable；
- integration 之后必须保留原始 assay / reduction；
- 需要同时报告 batch mixing 和 biological separation；
- AI 不能只根据 UMAP 外观判断 integration 成功。

### P1：marker、差异和 annotation

建议 action：

```text
find_markers
run_de_test
run_pseudobulk_de
score_gene_programs
map_cell_types
compare_annotations
```

关键约束：

- marker 结果必须保存 contrast、assay、layer、group definition；
- annotation 必须保留 reference、版本、方法和置信度；
- AI 生成的细胞类型名称不能覆盖原始 cluster identity；
- 低置信度 annotation 必须标记为候选而不是事实。

### P2：稀有细胞和异常群体

建议 action：

```text
inspect_rare_clusters          # summarize_small_cluster_evidence()（只聚合证据，不下结论）
run_doublet_diagnostic         # run_doublet_detection
run_rare_cell_detection        # run_rare_cell_detection（已实现）
compare_rare_cell_evidence     # 尚未实现
```

关键约束：

- 不允许仅按 cluster size 自动删除小群体；
- 需要同时检查 QC、doublet、sample replication 和 marker；
- 稀有群体必须报告“支持它的证据”和“替代解释”。

当前落地方式：`summarize_small_cluster_evidence()` 只聚合已经算好的信号并返回
`n_independent_signals_available`，自身不判断“真实/噪声”；`run_rare_cell_detection`
只写 `rare_cluster` 标签列，不删除、不过滤、不合并任何细胞；evidence 的 `claim_level`
按独立信号数量分级（0 个信号不登记 evidence，1 个信号 `associated` + `low_confidence = TRUE`，
≥2 个信号才 `consistent_with`）。

### P2：动态过程

建议 action：

```text
check_trajectory_inputs            # check_trajectory_readiness()（已实现）
run_trajectory                     # run_trajectory（已实现，start_cluster 必填）
compare_trajectory_roots           # 尚未实现
run_velocity                       # 尚未实现（需 spliced/unspliced）
compare_velocity_with_trajectory   # 尚未实现
run_fate_analysis                  # 尚未实现（依赖 velocity 输出）
```

关键约束：

- root、terminal state、direction 等关键假设必须显式记录；
- trajectory 不能把 cluster 顺序自动解释成时间顺序；
- velocity 不能在缺少 spliced / unspliced 的对象上伪造结果。

当前落地方式：`check_trajectory_readiness()` 只回答“有没有 cluster 和可用 embedding”，
**明确不给出任何起点建议**（“哪个 cluster 是起点”属于用户的生物学假设）；只读的
`summarize_trajectory_cluster_order()` 只描述 cluster 在嵌入维度上的分布，
`root_suggested` 恒为 `FALSE`。`run_trajectory` 包装 `RunSlingshot()`（而非
`RunSlingshot_trajectory()`），`group` 与 `start_cluster` 均为必填，
`prerequisites` 对缺失 / `NULL` / 空 / 不存在的 `start_cluster` 分别给出
`start_cluster_missing` / `start_cluster_unknown`，因此 AI 无法借 slingshot 的自动选根绕过确认。
evidence 以 `claim_level = "consistent_with"` 登记，只含 lineage 数量、pseudotime 分位数摘要和
匿名化后的 `cluster_N` 起点编码，不含逐细胞 pseudotime，并显式标注
`pseudotime_is_absolute_time = FALSE`。**本批未实现 velocity / CellRank / fate / spatial /
multimodal 的任何执行 action。**

### P3：程序、通路、通信和空间

建议 action：

```text
run_pathway_scoring
run_regulon_scoring
run_cell_cell_communication
run_spatial_niche
run_spatial_colocalization
run_multimodal_alignment
```

这些 action 需要更严格的依赖、版本、数据库和外部包环境管理，放在后续阶段。

---

## 8. 个性化分析计划 contract

高级 plan 不能只有一个动作列表，还需要表达“为什么选这条路线”和“如何判断路线是否成功”。

建议扩展 plan：

```r
list(
    plan_id = "integration_route_1",
    task = "batch_effect_diagnosis",
    user_goal = "保留 condition，降低技术 batch 影响",
    context_fingerprint = "...",
    assumptions = list(
        batch_is_technical = TRUE,
        condition_must_be_preserved = TRUE
    ),
    hypotheses = list(
        list(
            id = "batch_effect",
            statement = "PCA structure is partly explained by batch",
            evidence_required = c("pca_metadata_association", "sample_composition")
        )
    ),
    candidate_routes = list(
        list(id = "harmony", cost = "low", risk = "embedding_only"),
        list(id = "fastmnn", cost = "medium", risk = "parameter_sensitive"),
        list(id = "scvi", cost = "high", risk = "python_environment")
    ),
    selected_route = "harmony",
    success_criteria = list(
        "batch mixing improves",
        "condition separation is not erased",
        "marker structure remains interpretable"
    ),
    actions = list(...),
    requires_confirmation = TRUE
)
```

### 8.1 路线选择原则

AI 选择路线时必须同时考虑：

- 当前对象具备的输入；
- 用户保留或排除的变量；
- 统计和生物学目标；
- 计算成本；
- 外部依赖；
- 结果可解释性；
- 失败后的可恢复性；
- 是否存在可比较的替代路线。

### 8.2 不确定性要求

每个高级 recommendation 至少报告：

```r
list(
    recommendation = "...",
    confidence = "low|medium|high",
    evidence_refs = c(...),
    missing_evidence = c(...),
    alternative_explanations = c(...),
    user_decision_required = c(...)
)
```

---

## 9. Evidence graph 与 ledger 扩展

### 9.1 为什么只记录 state 不够

高级分析常常存在多条路线：

```text
raw → PCA → clustering
raw → Harmony → clustering
raw → fastMNN → clustering
raw → scVI → clustering
```

用户需要知道的不只是“哪个结果存在”，还包括：

- 哪些结果来自同一个 raw input；
- 哪些结果使用了不同校正方法；
- 哪些结论在多个路线中一致；
- 哪些结论只在一条路线出现；
- 哪些结论缺少样本层面的复现。

### 9.2 建议的 evidence node

```r
list(
    id = "evidence_001",
    kind = "deterministic_summary|plot|test|state|user_decision|ai_claim",
    source = "run_batch_diagnostic_1",
    claim_level = "observed|measured|associated|consistent_with",
    scope = list(
        cells = ..., 
        samples = ..., 
        groups = ...
    ),
    values = list(...),
    uncertainty = list(...),
    created_at = ...
)
```

### 9.3 AI claim 与 evidence 的关系

AI 不直接写入“事实”。它写入 claim：

```r
list(
    statement = "Cluster 4 is consistent with a cytotoxic lymphocyte population",
    claim_level = "consistent_with",
    evidence_refs = c("marker_result_4", "cluster_4_sample_distribution"),
    alternative_explanations = c("doublet", "ambient RNA"),
    confidence = "medium"
)
```

### 9.4 Evidence 引用验证和 claim ceiling

每个 `evidence_ref` 必须在当前 SCE 的 ledger/evidence registry 中解析到唯一记录，并验证：

- 来源 analysis/state 存在且状态为 completed；
- 证据对象、数据范围和 context fingerprint 与本次 claim 相符；
- 引用没有指向过期或已被替代的分析；
- claim 的 scope 不超出证据覆盖的 samples/groups/cells；
- AI claim 不可被当作独立 deterministic evidence 再次引用。

第一版允许的 claim level 为：

```text
observed | measured | associated | consistent_with
```

`causal` 不属于 AI 可生成的 claim level。只有用户明确提供、可追溯的外部因果研究设计和证据，并经过独立的因果分析流程后，才允许在报告层讨论因果解释；该解释仍不得由当前 AI action 自动生成或登记为 measured fact。

置信度 `low|medium|high` 是对证据充分性的定性摘要，不得冒充校准概率。若未来要输出概率，必须先定义校准数据集、校准方法和适用范围。

### 9.5 证据独立性

“独立证据”按数据来源和计算依赖定义，不按输出表的数量定义。若两个结果共享相同原始输入、相同聚类分组或相同 reference/database，则不能仅因为来自不同函数就视为独立。Evidence node 应带 `parents` 或 `dependency_group`；独立性判定由 R 根据 lineage 检查，AI 只解释该判定。


---

## 10. 人机协作边界

### 10.1 AI 可以自动做

- 读取和整理 ledger；
- 调用只读 diagnostics；
- 提出候选路线；
- 检查前置条件；
- 生成 dry-run；
- 解释 R 已经计算出的结果；
- 标记证据缺失和替代解释。

### 10.2 必须请求用户确认

- 删除或过滤细胞；
- 改变核心 biological variable；
- 选择 integration 变量；
- 选择 reference 和 annotation label；
- 选择 trajectory root 或 terminal state；
- 执行高成本或需要 Python 环境的分析；
- 覆盖已有 assay / reduction；
- 将候选 annotation 写入用户可见字段。

### 10.3 永远不能自动做

- 执行任意字符串形式的 R 代码；
- 删除原始 counts；
- 把 AI 推测写成 measured result；
- 没有证据时宣称因果关系；
- 在未知 batch / condition 语义时自动进行 integration；
- 在失败后把部分结果伪装为完成。

### 10.4 执行状态、部分失败和回滚语义

高级 plan 需要使用统一状态机：

```text
proposed → validated → dry_run → awaiting_confirmation
        → running → completed
        → failed | completed_with_errors | cancelled | stale
```

状态约束：

- `proposed` 不代表任何 action 已执行；
- `validated` 只代表当前 object、registry、参数和 prerequisite 通过检查；
- `dry_run` 不得改变 SCE；
- `awaiting_confirmation` 必须绑定 plan fingerprint、object fingerprint、action 参数和用户确认会话；
- `stale` 表示对象、plan、registry 或环境发生变化，必须重新验证；
- `completed_with_errors` 必须列出成功、失败和跳过的 action；
- 失败 action 不得将非法或不满足 contract 的返回值写入 `current`。

第一阶段不强制实现完整事务回滚，但必须准确区分：

1. **失败隔离**：失败返回值不会替换当前对象；
2. **部分完成**：前序成功 action 仍然存在，并在 ledger 中明确记录；
3. **事务回滚**：后续版本可通过对象快照或分支对象恢复全部执行前状态。

在完整回滚实现前，成功指标不得写成“失败不会污染对象”，应写成“失败不会写入非法返回值，并准确记录部分完成状态”。

### 10.5 设计语义 clarification

如果高级任务需要 `batch`、`condition`、`subject`、`sample`、`root`、`reference` 或 protected variables，而 profile 中无法确定其语义，AI 必须返回：

```r
list(
    status = "clarification_required",
    questions = list(...),
    assumptions = list(...),
    blocked_actions = list(...)
)
```

在用户回答前：

- 不得把候选 metadata 列直接当作事实；
- 不得执行依赖该语义的 action；
- 不得生成声称已完成的解释；
- 用户回答要作为 `user_decision` evidence 写回 ledger。

**当前实现状态（呈现层 + 记录层）**：判定逻辑本身未改动，本节只解决"如何呈现"和"如何留痕"。

- 呈现：`sclet_ai_format_clarification()` 把 `ValidateAIPlan()` 已经产生的错误字符串映射成结构化
  question 列表（design confirmation / annotation reference / annotation labels / trajectory root /
  reduction / counts assay），返回 `status` + `questions` + 原样保留的 `raw_errors`。它不做任何判定、
  不猜测答案、不推荐默认列；无法识别的错误按原文保留为 `unclassified` question，绝不丢弃。
  `RunAIAnalysis()` 在 `validation$valid == FALSE` 时把该结构挂到 `report$clarification`，
  `status` 仍为 `"invalid_plan"`（未改变取值集合）。
- 留痕：`sclet_ai_record_clarification_response()` 在用户回答后登记一条
  `kind = "user_decision"` 的 evidence，形成"AI 被阻塞 → 用户被提问 → 用户已回答"的审计链路。
  它**不**在 `ConfirmAIDesignSemantics()` 内部自动调用（确认机制本身保持原样），由上层 UX 流程显式串联。
  为遵守既有 `sclet_ai_evidence_value_ok()` 规则且不为其开特例，answer 的自由文本**不写入** evidence，
  只记录问题码、被点名 colData 列的序号以及若干布尔标记。
- 交互：`ResolveAIClarifications()` 把上述 question 逐条呈现给用户、读取答案、按需应用并留痕，
  至此"AI 被阻塞 → 用户被提问 → 用户已回答 → 可重新执行"形成闭环。它不代替用户作答：
  不推荐任何默认列 / cluster / reference；空答案记为 `skipped` 而不是被自动填上；
  指向不存在 colData 列的答案记为 `invalid`，既不写 confirmation 也不写 decision。
  只有 `design_batch` 会自动应用（经 `ConfirmAIDesignSemantics()`），其余问题描述的是 plan 参数，
  仅留痕并交还给调用方回填 plan（可用 `apply_answer` 钩子接管）。
  非交互会话且未提供 `ask` 时拒绝提示，直接返回 `needs_interactive` 且不改动对象。


## 11. API 演进建议

### 11.1 保留现有新手 API

```r
RunAIAnalysis(object, goal, confirm = "ask", ...)
AskAI(object, question, ...)
```

它们负责低门槛入口。

### 11.2 新增 profile API

```r
GetAIProfile(object, design = NULL, ...)
```

### 11.3 新增高级调查 API

```r
AIInvestigate(
    object,
    question,
    design = NULL,
    scope = c("diagnosis", "advanced", "report"),
    ...
)
```

`AIInvestigate()` 默认只读，输出 findings、evidence_refs、missing_evidence 和 candidate routes，不执行修改对象的 action。

### 11.4 新增多路线比较 API

```r
CompareAIAnalyses(
    object,
    ids = NULL,
    criterion = c("evidence", "stability", "biological_preservation", "cost"),
    ...
)
```

### 11.5 扩展 beginner facade

`RunAIAnalysis()` 后续增加：

```r
RunAIAnalysis(
    object,
    goal,
    design = NULL,
    mode = c("beginner", "advanced"),
    constraints = NULL,
    confirm = "ask",
    ...
)
```

其中：

- `mode = "beginner"`：隐藏复杂 plan 细节；
- `mode = "advanced"`：返回候选路线、证据、假设和完整 validation；
- `design`：提供 sample、batch、condition、subject 语义；
- `constraints`：提供成本、隐私、环境和生物学保护条件。

### 11.6 默认模型

所有上层 API 的 `model` 都保持可选：

```r
RunAIAnalysis(sce, goal = "...")
AskAI(sce, "...")
AIInvestigate(sce, "...")
```

解析优先级以当前实现为准：

1. 用户显式传入的 `model`（覆盖默认值）；
2. aisdk 已配置的默认模型；
3. 环境变量 `OPENAI_MODEL`；
4. 都未配置时返回清晰的 `sclet_ai_missing_model` 错误。

后续如需调整优先级，必须作为明确的 API 行为变更记录，而不能在不同入口间各自实现。

## 12. 高级分析成功标准与指标

高级 action 必须在 plan 中声明 `success_criteria`，并由 R 计算指标；AI 只能解释指标，不得自行判断“看起来更好”。第一版要求同时保留原始路线作为 baseline。

### 12.1 integration

至少报告：

- batch mixing：例如 graph LISI、kBET 或其他已选定的 batch mixing 指标；
- biological preservation：condition/label separation、marker/program preservation 或已知参考一致性；
- cluster stability：重复抽样或参数扰动下的 ARI/NMI、cluster overlap 或 marker overlap；
- cell/sample representation：是否删除或压低某个样本、condition 或稀有群体；
- runtime、内存和外部环境信息。

不允许把单一综合分数当作“integration 成功”。如果 batch mixing 改善但 biological preservation 下降，结果必须标记为 trade-off，并要求用户决定。

### 12.2 annotation

至少报告：

- reference、版本和映射方法；
- label-level confidence；
- marker/reference evidence；
- 未映射或冲突的细胞比例；
- 不同 reference 或方法之间的一致性；
- cluster identity 与 annotation label 的关系。

### 12.3 rare-cell 和异常群体

至少报告：

- 群体大小及其样本/subject replication；
- QC、doublet 和 ambient RNA 风险；
- marker 或 gene-program evidence；
- 是否只在一个样本或一个批次出现；
- 删除、合并和保留三种解释的比较。

### 12.4 指标 contract

每个指标记录：

```r
list(
    name = "graph_lisi",
    value = ..., 
    direction = "higher_is_better|lower_is_better|target_range",
    scope = list(samples = ..., groups = ..., cells = ...),
    method = ..., 
    parameters = ..., 
    baseline_ref = ..., 
    uncertainty = ..., 
    evidence_id = ...
)
```

如果指标不能计算，必须输出 `status = "not_available"` 和原因；不能用自然语言猜测替代。

---

## 13. 开发路线图

### Phase A：Dataset profiling 与诊断基础

**目标**：AI 能回答“这个对象是什么、缺什么、最值得检查什么”。

交付：

- `GetAIProfile()`；
- design 字段识别和用户确认；
- sample / batch / condition 摘要；
- QC by group；
- PCA / cluster 与 metadata 关联摘要；
- read-only diagnostic registry；
- profile schema 和测试。

完成标准：

- 不发送完整表达矩阵；
- 相同 SCE 得到稳定 profile fingerprint；
- profile 能明确区分 observed、missing 和 inferred；
- 至少覆盖 batch、rare cluster、annotation readiness 三类诊断。

### Phase B：高级分析路线规划

**目标**：AI 能基于 profile 选择分析路线，而不是只生成标准 pipeline。

交付：

- `AIInvestigate()`；
- hypothesis / evidence / alternative route contract；
- candidate route ranking；
- user decision checkpoint；
- advanced plan schema；
- 低成本 route preview。

完成标准：

- AI 不能在缺少 design 时自动猜 batch 和 condition；
- 每个 recommendation 都有 evidence_refs 和 missing_evidence；
- plan 能表达 success criteria 和 alternative routes。

### Phase C：批次、整合和稳定性

**目标**：优先解决最常见、最需要个性化判断的高级问题。

交付：

- batch diagnostics；
- integration action adapters；
- pre/post correction comparison；
- cluster stability；
- biological preservation checks；
- 多路线比较。

完成标准：

- 原始 assay 和 reduction 不被覆盖；
- 默认禁止校正 protected biological variables；
- 报告 batch mixing 与 biological separation 的共同变化；
- integration 结果可回溯到输入和参数。

### Phase D：annotation、DE 和稀有群体

**目标**：让 AI 能围绕细胞身份和小群体提出证据链，而不是直接命名。

交付：

- marker / DE / pseudobulk action；
- reference mapping；
- annotation confidence；
- rare-cell / doublet evidence chain（已实现第一批：`check_rare_cell_readiness()`、`summarize_small_cluster_evidence()`、`rare_cell` group 与按信号数量分级的 `claim_level`；decontX 独立信号与多路线比较 `compare_rare_cell_evidence` 仍待建设）；
- candidate label 与 confirmed label 的区分。

完成标准：

- AI 不能把候选 label 覆盖 cluster identity；
- annotation 必须记录 reference、方法和版本；
- 稀有群体必须有至少两类独立证据或明确标注低置信度。

### Phase E：动态、空间和多模态

**目标**：把 ledger + AI 扩展到复杂主线。

交付：

- trajectory / velocity / fate；
- pathway / regulon；
- communication / spatial niche；
- ADT / ATAC / multimodal；
- 跨主线 evidence graph。

完成标准：

- root、terminal state、reference 和 database 版本均有记录；
- 外部环境失败有 typed error；
- 高级结果可以被 AI 解释但不会被夸大为因果结论。

### Phase F：用户体验与评估

**目标**：把高级能力变成可用的研究助手，而不是一组底层函数。

交付：

- beginner / advanced 双模式；
- 对话式 QMD 和 vignette；
- 任务级回归数据集；
- plan validation benchmark；
- AI claim calibration benchmark；
- 成本、延迟、失败恢复报告。

---

## 14. 测试策略

### 14.1 确定性单元测试

覆盖：

- privacy policy 和 outbound payload redaction；
- profile schema；
- 设计字段识别；
- diagnostics 输出；
- action prerequisite；
- output/state contract；
- evidence refs；
- protected variables；
- 多路线比较；
- fingerprint 和 stale plan；
- failure isolation；
- `RunAIAnalysis()` beginner facade；
- 默认模型解析。

### 14.2 fixture 数据

至少准备以下小型 fixture：

1. 单样本 PBMC-like 数据；
2. 两批次、两 condition 数据；
3. 极小 rare population 数据；
4. 具有 spliced / unspliced 的 velocity fixture；
5. 带参考标签的 annotation fixture；
6. 含错误 metadata 和缺失 assay 的 negative fixture。

### 14.3 AI contract tests

使用 mock adapter 验证：

- AI 只看到 bounded context；
- outbound payload 不含完整矩阵、原始 ID、token 和未授权 metadata；
- 结构化结果可以规范化；
- evidence ref 必须解析到当前 ledger 中的 completed evidence；
- claim scope 不得超出 evidence scope；
- 缺失 evidence 时 claim level 不会升级；
- 非法 action 被拒绝；
- 缺少 design 时会返回 clarification request；
- stale plan、过期 confirmation 和 registry 变化会被拒绝；
- 部分失败会准确返回 completed_with_errors 或 failed；
- AI 不会把未执行的 action 写成 completed。

### 14.4 在线 smoke tests

继续采用双重门控：

```r
nzchar(Sys.getenv("DEEPSEEK_API_KEY"))
identical(Sys.getenv("SCLET_RUN_ONLINE_TESTS"), "true")
```

在线测试只验证：

- provider 连通；
- structured result 可解析；
- plan 输出能被 validation 处理；
- answer 中保留 uncertainty。

不把在线模型输出作为唯一回归标准。

---

## 15. 成功指标

### 用户体验

- 新手可以仅使用 `RunAIAnalysis()` 和 `AskAI()` 完成一次 guided analysis；
- 专家可以访问完整 plan、evidence 和 execution record；
- 用户不需要重复输入模型名称；
- AI 对象状态解释能够指出缺失前置条件。

### 科学可靠性

- AI 生成的每条高级 claim 都有 evidence_refs；
- 不允许 observed → causal 的无证据升级；
- 不允许候选 annotation 覆盖原始 identity；
- 多路线结果可以比较一致性和差异。

### 工程可靠性

- plan fingerprint 可验证；
- action 无法越过 registry；
- 失败不会写入非法返回值，并准确记录部分完成状态；完整事务回滚作为后续能力；
- 所有高级结果写入 lineage 和 ledger；
- 无 API key 时基础分析、离线测试和文档渲染正常。

### 性能与成本

- profile 不复制完整矩阵；
- 低成本 diagnostics 可以在本地快速完成；
- 高成本 action 在执行前显示预估成本；
- AI 不会重复执行已经完成且 fingerprint 未变化的分析。

---

## 16. 当前实现与本 Spec 的差距

当前已经完成：

- bounded ledger；
- structured AI result；
- read-only tools；
- native basic action catalog；
- plan validation；
- dry-run / confirmation；
- execution contracts；
- beginner `RunAIAnalysis()`；
- read-only `AskAI()`；
- 默认模型解析，不要求上层显式传 `model`。

还需要建设：

- 更完整的 sample / batch / condition 语义确认交互式 UX（`sclet_ai_format_clarification()` 结构化呈现、`sclet_ai_record_clarification_response()` 以 `user_decision` 留痕、`ResolveAIClarifications()` 交互式问答闭环均已实现；尚未实现的只是非终端形态的 UI，例如 Shiny / 网页端确认界面）；
- rare-cell / doublet 剩余项（第一批已实现，见上）：decontX ambient RNA 作为独立信号、`compare_rare_cell_evidence` 多路线比较、以及任何删除/合并稀有群体的 action（明确不在当前范围）；
- trajectory 剩余项（readiness + root-confirmed 执行已实现，见上）：`compare_trajectory_roots` 多 root 比较、`check_velocity_readiness`、`run_velocity`（需 spliced/unspliced）、`run_fate_analysis`（依赖 velocity 输出）；
- spatial 和 multimodal action catalog；
- hypothesis、success criteria 的 plan-level 执行与验证；
- advanced mode 的用户体验；
- 完整 causal claim ceiling（第一层已实现，见上：跨链 support 分组 + 确定性 ceiling）；
- causal / 机制性推论本身（目前 ceiling 只约束"证据支持强度"，不对"因果断言"设限，
  跨证据独立性审计也尚未接入 plan 或 result 层，仅作为独立只读 API 存在）；
- `run_integration` design_confirmed gate 方案 B（已实现）：独立 state `ai_design_confirmation` + 导出 API `ConfirmAIDesignSemantics(object, design = list(batch = "xxx", [condition = "...", subject = "..."]))`，删除 input_schema 里可伪造的 `.design_confirmed` 布尔；`prerequisites` 同时匹配三重条件才接受确认：① 存的骨架指纹（n_cells/n_features/colData 列名集合/是否有 rowData）与当前对象一致——仅结构变化才触发该层过期；② 存的 `summary$design_value_key`（design 中每个被点名 colData 列的逐单元值签名）与当前该列值逐格一致——同列名、但列下标签/值被重新赋义（如 batch 从 a/b 重映射为 T_cell/B_cell）时该层失效；③ `summary$design$batch == params$batch`。确认写完后再跑 PCA/写 preprocess/integration state 不会让确认失效（避免「确认一次就废」假阳性），但真正能影响 integration 语义的变更（batch 值的重赋义、列重命名/增删、对象维度变化）会让确认过期；从根源上杜绝 AI 自签通行证。

因此下一轮开发不应继续以“增加更多基础 action”为主，而应优先实现：

```text
GetAIProfile()
    ↓
高级 deterministic diagnostics
    ↓
AIInvestigate()
    ↓
批次 / integration / annotation / rare-cell 高级 action
    ↓
多路线比较和 evidence graph
    ↓
marker / DE / annotation 端到端证据链
```

---

## 17. 第一批建议立即开发的任务

第一批任务均已完成（含只读骨架 + 真实执行 + 隐私门控）：

1. **`GetAIProfile()`**：统一数据结构、实验设计候选字段、QC 和当前状态。
2. **确定性 diagnostics**：覆盖 QC 分组摘要、PCA–metadata 关联、cluster/sample 组成、小 cluster 和 integration readiness。
3. **`AIInvestigate()`**：只读地组合 profile、ledger 和 diagnostics；在设计语义不明确时返回 clarification requirement。

下一批（第 4–6 项）高级执行 Phase C 亦已完成：

4. **evidence graph integration**：`dependency_group` hash、`sclet_ai_evidence_independence()` lineage/dependency 查询、scope 越界与 cross-kind 互引拒绝。
5. **strict provider-side privacy**：`sclet_ai_call()` 默认 `enforce_privacy = TRUE`；deny 类永久硬拒绝不受 consent 影响；requires_user_consent 触发显式 `warning()` UX；metadata/payload 写入 `payload_policy`。
6. **真实 integration route execution**：`RunIntegrationRoutes()` 执行 raw baseline 登记 + fastMNN/Harmony/scVI；每条路线记录 batch mixing、biological preservation、cluster stability（`not_available` 带 reason）与 runtime；`CompareAIAnalyses()` 按 criterion 排序、detect batch-vs-bio tradeoff、缺失指标归一化为 `not_available`；recommendation 永远 NULL，不自动选优。

Phase C 收尾 + Phase D 第一批（第 7–8 项）已完成：

7. **design confirmation 硬化**：`ConfirmAIDesignSemantics()` 取代可被 AI 自称的 `.design_confirmed` 布尔；`run_integration` 的 `prerequisites` 按骨架指纹 + 逐列值签名双重校验确认记录，任一维度漂移（列增删改名、或同列名下取值被重新赋义）都会让旧确认过期。
8. **marker / DE / annotation 证据链**：`check_annotation_readiness()` 只读诊断；`AIDefaultExecutionRegistry` 新增 `annotation` group（`run_de_test`、`run_annotation`）；`run_annotation` 的 `ref`/`labels` 为必填参数，禁止落入 `RunSingleR()` 默认下载人类 reference 的分支；两个 action 的执行结果都登记为 evidence，`claim_level` 限定为 `associated`/`consistent_with`，不含逐细胞原始标签或基因名，且从不覆盖 `Idents()`。

此外已提供 `RunBasicWorkflow()` 基础 facade、结构化 clarification UX（`sclet_ai_format_clarification()` 呈现、`sclet_ai_record_clarification_response()` 以 `user_decision` evidence 留痕、`ResolveAIClarifications()` 交互闭环），以及跨证据 claim ceiling（`sclet_ai_claim_ceiling()`）。

`claim ceiling` 的设计要点：evidence 节点先按"独立支撑线"分组，再算 ceiling。两个节点若共享
`dependency_group`、存在 `parents` 祖先关系、或来自同一次底层分析（`source` 相同），则视为同一
条支撑线；再取该关系的连通分量作为分组，保证独立性只会被低估、不会被高估。这一条在实践中很关键：
单次 rare-cell 运行会为每个小群体各产生一个 evidence 节点，按节点计数会把"一次测量"夸大成
"几十个独立信号"。ceiling 规则是确定性的：无支撑则不可断言；单条支撑线封顶 `associated`，
孤立节点无法被升级；≥2 条独立支撑线时 ceiling 介于 `associated` 与 `consistent_with` 之间，
并受最弱的那个节点约束。`user_decision` 节点不计入支撑数，与既有"AI 生成的 claim 不得引用
人类决策"规则一致。该函数只做评估，不生成、不强化、不登记任何 claim。

当前版本已提供 evidence registration、payload sanitizer、真实多路线执行、真实路线比较、默认 privacy gate、lineage/dependency 独立性查询、design confirmation 硬化、marker/DE/annotation 证据链、rare-cell/doublet 证据链（`check_rare_cell_readiness()`、`summarize_small_cluster_evidence()`、`rare_cell` group、按独立信号数量分级的 `claim_level`）与 trajectory 第一批（`check_trajectory_readiness()`、只读 `summarize_trajectory_cluster_order()`、`trajectory` group、root 必填的 `run_trajectory`）；下一阶段建议建设 velocity/CellRank/fate、spatial/multimodal action catalog、trajectory 多 root 比较和 rare-cell 多路线比较，使 sclet AI 从"integration + annotation 审计助手"进一步扩展为"全流程高级分析研究助手"。
