# Handoff：Marker / DE / Annotation 证据链（Phase D 第一批）

- **交接对象**：另一个 AI coding agent
- **仓库**：`/home/wang/data/source/omics/sclet`，分支 `devel`
- **前置任务状态**：Phase A（profile/diagnostics）、Phase B（AIInvestigate）、Phase C（evidence graph、privacy、integration route execution/comparison、design confirmation 硬化）均已验收通过。本任务在此基础上做**下一环**：把已有的 marker/DE/annotation 底层函数接成 AI-native 证据链。
- **不要动**：`R/ai-comparison.R`、`R/ai-evidence.R`、`R/ai-privacy.R`、`R/ai-integration-routes.R`、`R/ai-design-confirmation.R`、`run_integration` action。这些已验收，出问题会导致回归。

---

## 0. 一句话目标

现在 sclet AI 已经能真实执行 integration 路线、登记 evidence、比较路线。但 **annotation（细胞类型标注）和 marker/DE（差异表达）这条线还是空的**——底层算法早就存在（`RunSingleR()`、`FindMarkers()`、`RunDEtest()`、`FindAllMarkers()`），只是没有 AI-native 层：没有 action 注册、没有 evidence 登记、没有"候选标注 vs 确认标注"的区分。

你要做的是把这条线补上，**遵循 integration 那一套已经验证过的模式**：deterministic diagnostics → action registry（受控执行）→ evidence 登记 → 候选/确认区分。

---

## 1. 硬约束（和之前几轮一样，不要违反）

```bash
# 日常验证
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'

# 完整检查，必须用 Makefile
make check
```

- `make check` 必须 `0 errors | 0 warnings | 0 notes`；
- 不要直接跑 `R CMD check`；
- 不要跑 `make rd`（会重写 NAMESPACE，丢手工 import）；手工维护 `NAMESPACE`、`DESCRIPTION` 的 `Collate:`、`man/*.Rd`；
- `R/*.R` 必须纯 ASCII（中文用 `\uXXXX`，或者更简单：**新代码里的注释和字符串直接写英文**，不要冒非 ASCII 的风险）；
- 不要泄露密钥；
- 不要破坏现有 335+ 个通过的 AI 测试；
- 提交前跑：

```bash
cd /home/wang/data/source/omics/sclet
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"
git diff --check
```

**本机依赖现状**（实测）：`SingleR`、`celldex`、`batchelor` 都已安装，`harmony` 未安装。因此 annotation 相关测试可以真实跑通 `RunSingleR()`（不需要 mock），但如果你的测试涉及需要下载参考数据集的场景（比如 `celldex::HumanPrimaryCellAtlasData()` 联网拉取），必须 mock 掉，不要让测试依赖网络。

---

## 2. 现有底层函数（不要重新实现，只包一层 AI-native adapter）

| 函数 | 文件 | 作用 | 已有 state type |
|---|---|---|---|
| `RunSingleR(object, ref=NULL, labels=NULL, ...)` | `R/annotation.R:15` | reference-based 细胞类型标注，写入 `colData` 的 `*_labels`/`*_pruned.labels`/`*_score` 列 | `annotation` |
| `RunReferenceMapping()` | `R/annotation.R:160` | 另一种 reference mapping | `annotation` |
| `RunKNNPredict()` | `R/annotation.R:221` | KNN 标签迁移 | `annotation` |
| `RunSymphonyMapping()` | `R/annotation.R:372` | Symphony reference mapping | `annotation` |
| `FindMarkers(object, ident.1=NULL, ident.2=NULL, ...)` | `R/markers.R:142` | 两组之间的差异表达 | 无独立 state，走 `RunDEtest` |
| `FindAllMarkers()` | `R/markers.R:21` | 每个 cluster 的 marker | 同上 |
| `RunDEtest(object, ident.1=NULL, ident.2=NULL, name="detest", ...)` | `R/markers.R:178` | 统一入口，登记 `detest` state，`artifacts$result` 存结果 data.frame | `detest` |

`sclet_state_types()`（`R/state.R`）已经包含 `"annotation"` 和 `"detest"`，**不需要新增 state type**。

---

## 3. 要实现的 4 件事

### 3.1 确定性 annotation/marker readiness 诊断

在 `R/ai-diagnostics.R` 里新增（跟 `check_integration_readiness()` 同一个文件、同一种风格）：

```r
check_annotation_readiness(object, design = NULL)
```

要求：

- 检查对象是否有可用的 assay（跟 `check_integration_readiness` 一样）；
- 检查是否已有 cluster/ident 信息（`Idents(object)` 或类似）——marker/DE 需要分组，没有分组就不能跑；
- **不要**自动猜测哪个 reference 数据库适用于这个物种/组织——如果 `design` 或某种显式参数里没有指明 reference 来源，返回 `clarification_required`，列出需要用户回答的问题（例如"这是人类数据还是小鼠数据？""希望用哪个 reference？"）；
- 返回结构参考 `check_integration_readiness()` 的风格：`status`（`ready_for_diagnostic` / `clarification_required` / `not_ready`）、`checks`、`blocked_actions`、`questions`。

同时把 `R/ai-profile.R` 里现有的 `annotation_readiness`（纯粹按列名猜的 stub，第 416 行）**保留不动**——那是 profile 里的轻量候选提示，不要跟这个新的确定性诊断函数混在一起，两者服务的场景不同（profile 是"有哪些候选列名"，这个新函数是"能不能真的跑 annotation"）。

### 3.2 Action registry：新增 `annotation` group

在 `R/ai-execution.R` 的 `AIDefaultExecutionRegistry()` 里：

- `allowed_groups` 加入 `"annotation"`（参考第 178 行 `integration` 是怎么加进去的）；
- 新增至少两个 action：
  - `run_de_test`：包 `RunDEtest()`；
  - `run_annotation`：包 `RunSingleR()`（优先选它，因为本机已装依赖，测试能真实跑）。

**每个 action 都要有 `prerequisites`**，并遵循和 `run_integration` 完全一样的安全模式：

- `run_de_test` 的 `prerequisites`：
  - 检查 `ident.1`/`ident.2`（如果指定）对应的分组是否存在；
  - 如果 `ident.1 = NULL`（跑 `FindAllMarkers`），检查当前是否有 active identity/cluster；
  - **不需要** design confirmation gate（marker/DE 本身不像 integration 那样涉及"批次语义被误判"的风险，跑错了组别顶多是统计结果没意义，不会造成"AI 悄悄纠正生物学变量"这种更严重的问题）——但要检查最小样本量（每组至少若干个细胞，避免统计上没有意义的比较）。

- `run_annotation` 的 `prerequisites`：
  - 检查 `ref`/`labels` 参数是否提供（不能让 AI 自己决定用哪个 reference——这是本任务真正的风险点，处理方式见 3.3）；
  - 依赖检查：`requireNamespace("SingleR")`/`requireNamespace("celldex")`，缺失时返回可读原因（模仿 `run_integration` 对 `harmony` 缺失的处理，`R/ai-execution.R` 第 517-521 行附近）。

`estimated_cost`：`run_de_test` 用 `"medium"`（Presto 计算有一定开销），`run_annotation` 用 `"medium"`（SingleR 分类）。`requires_confirmation = TRUE`（两者都会修改对象的 colData/state）。`allowed_state_types`：`run_de_test` 是 `"detest"`，`run_annotation` 是 `"annotation"`。

默认 registry（`include = "read"`）不应该包含这两个 action，必须显式 `include = c("read", "annotation")` 才能用到。

### 3.3 Reference 来源必须显式，不能由 AI 猜

这是本任务**最容易做错的一步**，务必仔细看：

Spec 里已经写死的规则（`.dev/ai-advanced-analysis-spec.md` 第 12.2 节 annotation 部分）：annotation 必须记录 reference、方法、版本、置信度；AI 生成的细胞类型名称不能覆盖 cluster identity；低置信度必须标记为候选而非事实。

具体到 `run_annotation` action：

- `input_schema` 里 `ref` 和 `labels` 应该是**必填参数**（`required = TRUE`），或者提供一个受限的枚举（例如只允许 `"HumanPrimaryCellAtlasData"`、`"MouseRNAseqData"` 这类明确命名的 celldex 数据集标识符，由 handler 内部去 `celldex::` 里查找对应函数）；
- **不要**让 `ref = NULL` 走到 `RunSingleR()` 内部的默认值逻辑（`R/annotation.R` 第 21-27 行——那里默认会去下载 `HumanPrimaryCellAtlasData`，这是"没有明确说是人类数据，就默认当成人类数据"，这正是不该由 AI 自动决定的假设）。在 action 的 `prerequisites` 或 `handler` 里显式拒绝 `ref` 为空的调用。

### 3.4 Evidence 登记 + candidate/confirmed 标签区分

annotation 结果不能直接当作事实写回。执行完 `run_annotation` 之后：

- 调用 `RecordAIEvidence()`（`R/ai-evidence.R`，上一轮已实现）把标注结果登记为一条 evidence：
  - `kind = "deterministic_summary"`；
  - `source` 指向刚写入的 `annotation` state id；
  - `values` 里只放**聚合摘要**（比如每个 cluster 落进哪个候选标签的比例、平均 confidence score），**不要**把逐细胞的 label 明细塞进 evidence（那属于原始/半原始数据，参考上一轮 `sclet_ai_evidence_value_ok()` 的边界，逐细胞标签数组长度可能很大且接近"原始数据"）；
  - `claim_level = "consistent_with"`（不能是 `"measured"` 或 `"observed"`——细胞类型标注本质是模型推断，不是直接测量）。
- 在 action 的 `output_schema` 或返回结果里，明确区分：
  - `colData` 里写入的原始预测列（`*_labels`/`*_pruned.labels`）——这是 `RunSingleR()` 已有行为，保留；
  - 一个新的、AI 层面的"候选标注"概念——**不要**用标注结果覆盖或重命名任何已有的 cluster identity/Idents 列。如果你打算新增一个便于用户看的汇总（比如"cluster 3 大概率是 T 细胞，置信度 0.82，证据 id xxx"），把它放在 evidence 或 action 返回值里，不要直接改 `Idents(object)`。

### 3.5（可选，时间不够可以跳过）marker/DE 结果的证据登记

如果时间允许，`run_de_test` 执行完之后同样调用 `RecordAIEvidence()`，`values` 里放聚合摘要（比如显著 marker 数量、top N marker 基因名、检验方法），不要把完整的 marker data.frame 塞进去。这条不是本轮强制项，但如果做了要在报告里说明。

---

## 4. 必须补的测试

### 4.1 readiness 诊断

```r
test_that("check_annotation_readiness requires cluster assignment", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    result <- check_annotation_readiness(sce)
    expect_true(result$status %in% c("not_ready", "clarification_required"))
})
```

### 4.2 action registry 隔离

```r
test_that("AIDefaultExecutionRegistry has annotation group only when explicitly included", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    default_registry <- AIDefaultExecutionRegistry(sce)
    expect_false("run_annotation" %in% names(default_registry))
    full_registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    expect_true("run_annotation" %in% names(full_registry))
})
```

### 4.3 reference 不能为空——真正验证拒绝逻辑

```r
test_that("run_annotation rejects calls without an explicit reference", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(rpois(40,5), 10L, 4L)))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    action <- registry$run_annotation
    result <- action$prerequisites(sce, list())  # 没传 ref/labels
    expect_true(is.character(result))  # 应该返回可读拒绝原因，不是 TRUE
})
```

### 4.4 端到端：真正走 `ValidateAIPlan` + `ExecuteAIPlan`

参考上一轮 `run_integration` 的端到端测试写法（`tests/testthat/test-ai-execution.R` 里搜 `"standard plan pipeline"`），给 `run_de_test` 写一条类似的、真正经过标准管线的测试，不要只测 `prerequisites()`。`run_de_test` 不依赖 harmony/celldex，应该能不 mock 真实跑通（用小型模拟数据）。

### 4.5 Evidence 登记验证

```r
test_that("run_annotation results are recorded as consistent_with evidence, not measured fact", {
    # 执行 run_annotation（真实跑或按需 mock celldex 下载），
    # 断言产生的 evidence node 的 claim_level == "consistent_with"，
    # 且 evidence$values 里不包含逐细胞的完整标签数组（只有聚合摘要）
})
```

---

## 5. 验收标准

- §4 五条测试类别全部有对应测试且通过（4.5 如果 3.5 跳过，可以只测最小闭环）；
- `run_annotation` 在没有显式 `ref`/`labels` 时，通过 `ValidateAIPlan()` 或 `prerequisites()` 被拒绝，理由可读；
- 默认 registry 不包含 `run_de_test`/`run_annotation`；
- annotation 产生的 evidence `claim_level` 不是 `measured`/`observed`；
- 不修改任何已验收模块（见文档开头"不要动"列表）；
- 跑：

```bash
cd /home/wang/data/source/omics/sclet
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"
git diff --check
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'
make check
```

期望：

```text
[ FAIL 0 | WARN 0 | SKIP 1 | PASS >= 340 ]
0 errors | 0 warnings | 0 notes
```

---

## 6. 交付物

- [ ] `R/ai-diagnostics.R` 新增 `check_annotation_readiness()`
- [ ] `R/ai-execution.R` 新增 `annotation` group、`run_de_test`、`run_annotation` 两个 action
- [ ] 对应的手写 `man/*.Rd`（新增导出符号都要有文档）
- [ ] `NAMESPACE` 新增 export（不要有重复行）
- [ ] `DESCRIPTION` 的 `Collate:` 如果新增了文件要同步
- [ ] `NEWS.md` 顶部加一条本阶段条目
- [ ] `.dev/ai-advanced-analysis-spec.md` 更新：把"marker / DE / annotation evidence chain"从"未实现"移到"已实现"，并简述做了什么、还差什么
- [ ] §4 的测试全部添加并通过
- [ ] 最终报告：贴出 §5 四条命令的真实输出

## 7. 不要做的事

- 不要让 `ref = NULL` 时静默使用默认的人类 reference（`HumanPrimaryCellAtlasData`）——这正是本任务要堵的"AI 替用户做了一个未声明的生物学假设"的洞。
- 不要用标注结果直接覆盖 `Idents(object)` 或任何现有 cluster identity 列。
- 不要把逐细胞的完整标签数组塞进 evidence `values`。
- 不要给 annotation 结果打 `claim_level = "measured"` 或 `"observed"`。
- 不要动已验收的 integration/evidence/privacy/design-confirmation 模块。
- 不要跑 `make rd`。

## 8. 完成后

把最终报告（§5 四条命令的真实输出）留下。审核会重点检查：`ref=NULL` 时是否真的被拒绝（不是只在文档里说拒绝）、evidence 的 `claim_level` 是否真的不是 `measured`、默认 registry 是否真的不包含这两个新 action、以及是否有任何一条新测试是"绕开标准管线直接调内部函数"的假端到端测试（上两轮都出现过这个问题，这次要避免）。
