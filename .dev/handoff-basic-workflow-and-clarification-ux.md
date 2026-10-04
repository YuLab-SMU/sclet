# Handoff：`RunBasicWorkflow()` facade + 结构化 clarification UX

- **交接对象**：另一个 AI coding agent
- **仓库**：`/home/wang/data/source/omics/sclet`，分支 `devel`
- **前置任务状态**：Phase C/D 已实现的模块（integration route execution、design confirmation 硬化、marker/DE/annotation 证据链、rare-cell/doublet 证据链——如果这些还没全部并入你的工作分支，以 `.dev/ai-advanced-analysis-spec.md` 里的"已实现"清单为准）。
- **不要动**：`R/ai-comparison.R`、`R/ai-evidence.R`、`R/ai-privacy.R`、`R/ai-integration-routes.R`、`R/ai-design-confirmation.R`、任何 `run_*` action 的 `prerequisites` 逻辑本身、`check_*_readiness()` 的判定逻辑本身。**本任务不改变任何安全边界的判定条件，只改变"判定结果如何呈现给用户"这一层。**

---

## 0. 两个任务，为什么放一起

这是 Spec 里"还需要建设"清单的第 2、3 项，都属于**用户体验层**的工作，不涉及新的安全边界：

1. **`RunBasicWorkflow()`**：新手友好的基础流程 facade，纯确定性 R workflow，不需要 AI 参与决策。之前一直是"拟议接口"没人实现。
2. **结构化 clarification UX**：现在 `RunAIAnalysis()`/`AIInvestigate()` 遇到"design 语义未确认"这类情况时，只是把校验错误的原始字符串（比如 `"design_semantics_not_confirmed: call ConfirmAIDesignSemantics(...)"`）原样抛给用户，不够结构化、不够友好，而且**用户回答后的确认动作不会被记录成 evidence**——Spec 第 10.5 节写了"用户回答要作为 `user_decision` evidence 写回 ledger"，这一步现在完全没做。

两者都是纯粹的"包装/呈现"工作，风险很低，适合一起做。**关键约束是：不要重新实现或修改任何 action 的 `prerequisites`、任何 `check_*_readiness()` 的判定逻辑、任何 `ConfirmAIDesignSemantics()`/`RecordAIEvidence()` 的校验规则**——这些都已经过审核，本任务只是在它们之上加一层更好的呈现和记录。

---

## 1. 硬约束（跟之前几批一样）

```bash
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'
make check
```

- `make check` 必须 `0 errors | 0 warnings | 0 notes`；
- 不要直接跑 `R CMD check`；不要跑 `make rd`；
- `R/*.R` 必须纯 ASCII；
- 不要破坏现有测试（用 `devtools::test(filter="ai-")` 先跑一遍看当前基线数字，完成后数字应该只增不减）；
- 完成后清理你自己验证用的临时脚本文件，不要留在仓库里；
- 提交前跑：

```bash
cd /home/wang/data/source/omics/sclet
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"
git diff --check
```

---

## 2. 任务一：`RunBasicWorkflow()`

### 2.1 目标

一个纯确定性、不需要 AI 参与的新手 facade，把标准单细胞基础流程串起来：

```r
sce <- RunBasicWorkflow(sce)
```

对应的底层函数链（都已存在，直接调用，不要重新实现）：

| 步骤 | 函数 | 文件 |
|---|---|---|
| 1 | `NormalizeData(object, scale.factor=10000, assay="counts")` | `R/preprocessing.R:56` |
| 2 | `FindVariableFeatures(object, nfeatures=2000, method="scran", ...)` | `R/preprocessing.R:139` |
| 3 | `ScaleData(object, features=NULL, assay="logcounts")` | `R/preprocessing.R:693` |
| 4 | `RunPCA(object, subset_row=NULL, exprs_values=NULL, layer=NULL, ncomponents=50, ...)` | `R/dimred.R:233` |
| 5 | `FindNeighbors(object, dims, reduction=NULL, k=10)` | `R/clustering.R:11` |
| 6 | `FindClusters(object, resolution=0.5)` | `R/clustering.R:79` |
| 7 | `RunUMAP(object, dims=NULL, reduction=NULL, layer=NULL)` | `R/dimred.R:306` |

### 2.2 设计要求

```r
RunBasicWorkflow(
    object,
    n_features = 2000,
    n_pcs = 30,
    cluster_resolution = 0.5,
    steps = c("normalize", "variable_features", "scale", "pca", "neighbors", "clusters", "umap"),
    verbose = TRUE
)
```

- 所有参数都有合理默认值，**这是一个"不用管细节就能跑通"的入口**，不是给专家用的可配置管线（专家应该直接调底层函数，不需要 `RunBasicWorkflow()`）；
- `steps` 参数允许跳过某些步骤（比如已经跑过 QC 和 normalize，只想从 PCA 开始）——每个 step 名对应上表的一行，实现时用一个 named list 把 step 名映射到对应的函数调用；
- **不需要**做任何 AI 调用、不需要接入 `AIAction`/registry——这是纯 R 函数，不属于 AI-native 那一层；
- `verbose = TRUE` 时用 `message()` 报告当前在跑哪一步（不要用 `cat()`，参考包里其它函数的日志风格，检查 `R/preprocessing.R`/`R/dimred.R` 里已有函数是怎么打印进度的，保持风格一致）；
- 每一步底层函数已经自己处理状态登记（`NormalizeData`/`RunPCA`/`RunUMAP`/`FindNeighbors`/`FindClusters` 调用 `sclet_set_analysis_state()`；`FindVariableFeatures`/`ScaleData` 只调用较轻量的 `sclet_log_command()`，不写完整 analysis state——这是现有行为，不需要你去补，`RunBasicWorkflow()` 只需要顺序调用，不用重复登记任何东西）；
- 返回更新后的 `object`，不需要返回任何 plan/validation/execution 结构（这跟 `RunAIAnalysis()` 不是同一类接口，不要照抄它的返回结构）。

### 2.3 需要你自己确认的点

- `FindNeighbors(object, dims, ...)` 的 `dims` 参数没有默认值（必填），你需要在 `RunBasicWorkflow()` 里根据 `n_pcs` 构造 `dims = 1:n_pcs` 传进去，先确认这样传参是否符合 `FindNeighbors()` 的预期用法（读一下它的实现和现有调用它的地方，比如之前几批 handoff 里的测试代码是怎么调的）；
- 确认 `RunPCA()` 的 `ncomponents` 参数跟 `RunBasicWorkflow()` 的 `n_pcs` 是同一个含义，不要假设。

### 2.4 测试

```r
test_that("RunBasicWorkflow runs the standard pipeline end-to-end with defaults", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(50 * 40, lambda = 5), nrow = 50L, ncol = 40L))
    )
    result <- RunBasicWorkflow(sce, n_pcs = 5)
    expect_true("PCA" %in% SingleCellExperiment::reducedDimNames(result))
    expect_true("UMAP" %in% SingleCellExperiment::reducedDimNames(result))
    expect_false(is.null(ActiveIdent(result)))
})

test_that("RunBasicWorkflow respects a restricted steps argument", {
    # 只传 steps = c("normalize", "variable_features")，
    # 断言没有产生 PCA/UMAP reduction
})
```

---

## 3. 任务二：结构化 clarification UX

### 3.1 现状问题

现在遇到未确认语义的情况，用户看到的是这种体验：

```r
result <- RunAIAnalysis(sce, goal = "帮我做 batch integration")
# result$status == "invalid_plan"
# result$report$errors == c("prerequisites not met for action integrate: design_semantics_not_confirmed: call ConfirmAIDesignSemantics(...)")
```

用户拿到的是一句给开发者看的错误字符串，不是一个"这是需要你回答的问题清单"的结构。而且，就算用户去调用了 `ConfirmAIDesignSemantics()`，这次"人类确认了 batch 语义"的事件本身**不会被记录成一条 evidence**——它只是写进了 `ai_design_confirmation` state（这个是对的，不要动），但没有跟"AI 曾经因为这个原因被阻塞"这件事关联起来，形成一个可追溯的"AI 问了问题 → 用户回答了 → 继续执行"的链路。

### 3.2 要做的事：新增一层展示/记录函数，不改判定逻辑

**新增**（建议放在 `R/ai-functions.R`，跟 `RunAIAnalysis()`/`AskAI()` 放一起）：

```r
sclet_ai_format_clarification <- function(validation_errors, readiness_results = list())
```

这是一个纯格式化函数：接收 `ValidateAIPlan()` 返回的 `errors`（字符串向量）和可选的、已经调用过的 `check_*_readiness()` 结果，**从已有的错误信息里解析/映射**出结构化问题列表，不重新做任何判定。返回结构参考 Spec 第 10.5 节：

```r
list(
    status = "clarification_required",
    questions = list(
        list(
            id = "design_batch",
            text = "AI wanted to run integration but the technical batch column has not been confirmed. Which colData column represents the technical batch? Call ConfirmAIDesignSemantics(object, design = list(batch = '<column name>')) to confirm.",
            blocked_action = "run_integration",
            related_function = "ConfirmAIDesignSemantics"
        )
        # 每个被拒绝的 action 对应一条 question
    ),
    raw_errors = validation_errors  # 保留原始信息，不要丢弃，专家用户可能想看
)
```

具体怎么从错误字符串里识别出"这是一个 design confirmation 问题"还是"这是一个 reference 缺失问题"还是"这是一个 root 未指定问题"——可以用简单的关键词匹配（`grepl("design_semantics_not_confirmed", ...)`、`grepl("reference_missing", ...)`、`grepl("start_cluster_missing", ...)`），因为这些错误信息本身已经是结构化命名的（前缀是固定的），不需要复杂的 NLP。**如果关键词匹配不到已知模式，把它当作一条通用问题原样保留，不要丢弃或报错。**

### 3.3 接入 `RunAIAnalysis()`

`RunAIAnalysis()` 目前 `validation$valid == FALSE` 时的返回块（`R/ai-functions.R` 约第 232-251 行）里，把 `report` 字段改成调用 `sclet_ai_format_clarification()` 产出的结构，**不要改变 `status` 字段的值仍然是 `"invalid_plan"`**（保持向后兼容，已有测试断言这个值），只是让 `report` 里多一层结构化信息：

```r
report = list(
    status = "invalid_plan",
    goal = goal,
    errors = validation$errors,
    warnings = validation$warnings,
    clarification = sclet_ai_format_clarification(validation$errors)
)
```

### 3.4 用户确认后，把这次确认记录为 evidence

这是本任务最重要的一步。在 `ConfirmAIDesignSemantics()` 的调用点之外（**不要修改 `ConfirmAIDesignSemantics()` 本身**——它已经验收，只负责写 `ai_design_confirmation` state），新增一个可选的、显式调用的辅助函数：

```r
sclet_ai_record_clarification_response <- function(object, question_id, answer, blocked_action = NULL)
```

调用 `RecordAIEvidence()`，`kind = "user_decision"`（这是已经允许的 evidence kind，不需要改 `RecordAIEvidence()` 本身），`values` 里放 `question_id`、`answer`（如果 answer 是一个 colData 列名这种非敏感标识符可以直接记录；如果 answer 里可能包含更敏感的内容，按 `sclet_ai_evidence_value_ok()` 现有规则该拒绝的地方让它自然拒绝，不要在这里加特例绕过）。

这个函数**不是**在 `ConfirmAIDesignSemantics()` 内部自动调用的（design confirmation 那条线已经验收，不要在里面加新逻辑），而是给上层 UX 流程用的：比如未来一个交互式脚本可以先展示 `sclet_ai_format_clarification()` 的问题，用户回答后，脚本先调用 `ConfirmAIDesignSemantics()`（走既有确认机制），再调用这个新函数把"用户被问过什么、回答了什么"这件事额外记一条 evidence，方便后续审计"AI 曾经因为什么原因被阻塞过、用户是怎么回应的"。

**如果你觉得这两步应该合并成一个函数**（调用一次就同时做 `ConfirmAIDesignSemantics()` + 记录 evidence），可以提议这个设计，但要在报告里说明为什么，并确保**不修改 `ConfirmAIDesignSemantics()` 的现有签名和行为**（可以新增一个包装函数，不要改原函数）。

### 3.5 测试

```r
test_that("sclet_ai_format_clarification maps a design-confirmation error to a structured question", {
    errors <- c("prerequisites not met for action integrate: design_semantics_not_confirmed: call ConfirmAIDesignSemantics(object, design = list(batch = 'batch', ...)) before running integration")
    result <- sclet_ai_format_clarification(errors)
    expect_equal(result$status, "clarification_required")
    expect_true(length(result$questions) >= 1L)
    expect_equal(result$raw_errors, errors)
})

test_that("RunAIAnalysis surfaces a structured clarification report on invalid_plan", {
    # 构造一个会触发 design_semantics_not_confirmed 的场景，跑 RunAIAnalysis()
    # 断言 result$status == "invalid_plan"（不变）
    # 断言 result$report$clarification$status == "clarification_required"
    # 断言 result$report$clarification$questions 非空
})

test_that("sclet_ai_record_clarification_response writes a user_decision evidence node", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    sce <- sclet_ai_record_clarification_response(sce, question_id = "design_batch",
        answer = "batch", blocked_action = "run_integration")
    ledger <- GetAnalysisLedger(sce, detail = "summary", include_artifacts = FALSE, include_data = FALSE)
    ev <- ledger$state_records$ai_evidence
    found <- Filter(function(x) identical(x$summary$kind, "user_decision"), ev)
    expect_true(length(found) >= 1L)
})

test_that("clarification formatting does not alter ValidateAIPlan or prerequisites behavior", {
    # 回归性检查：构造跟审核 run_integration 时同样的攻击场景
    # （AI 自称 .design_confirmed 或省略必填参数），
    # 确认 ValidateAIPlan 的 valid/errors 跟本任务开始前完全一致，
    # 只是 RunAIAnalysis 的 report 层多了 clarification 字段
})
```

最后一条测试很重要：它是用来防止你在做"呈现层"改进的时候不小心动到了判定逻辑本身。

---

## 4. 验收标准

- `RunBasicWorkflow()` 能跑通标准流程，`steps` 参数能正确跳过步骤；
- `sclet_ai_format_clarification()` 能把已知的错误模式（design/reference/root 三类）映射成结构化问题，未知模式原样保留不丢失信息；
- `RunAIAnalysis()` 的 `status` 字段值保持不变（仍是 `"invalid_plan"`），只是 `report` 多了 `clarification` 字段；
- 新增 `sclet_ai_record_clarification_response()`，能把用户确认动作登记为 `kind="user_decision"` 的 evidence；
- **回归测试确认 `ValidateAIPlan()`/`prerequisites` 的判定结果没有被这批改动影响**（§3.5 最后一条测试）；
- 不修改 `ConfirmAIDesignSemantics()`、任何 `check_*_readiness()`、任何 action 的 `prerequisites` 本身；
- 跑：

```bash
cd /home/wang/data/source/omics/sclet
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"
git diff --check
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'
make check
```

期望：`FAIL 0 | WARN 0`，`make check` → `0 errors | 0 warnings | 0 notes`。

---

## 5. 交付物

- [ ] `RunBasicWorkflow()` 新函数（放在 `R/` 下合适的文件，比如新建 `R/basic-workflow.R`，同步更新 `DESCRIPTION` 的 `Collate:`）
- [ ] `sclet_ai_format_clarification()`、`sclet_ai_record_clarification_response()`（放在 `R/ai-functions.R`）
- [ ] `RunAIAnalysis()` 的 `report` 字段接入 clarification 结构
- [ ] 对应的手写 `man/RunBasicWorkflow.Rd`（这个函数要 export）
- [ ] `NAMESPACE` 新增 `export(RunBasicWorkflow)`（`sclet_ai_format_clarification`/`sclet_ai_record_clarification_response` 是否需要 export，参考包内命名习惯——加了 `sclet_` 前缀的通常是内部函数，不导出，但如果你觉得用户应该能直接调用第二个函数去手动记录确认历史，可以导出，在报告里说明理由）
- [ ] `NEWS.md` 顶部加一条本阶段条目
- [ ] `.dev/ai-advanced-analysis-spec.md` 更新：`RunBasicWorkflow()` 从"拟议接口，尚未实现"改为已实现；第 10.5 节补充说明 clarification 现在有了结构化呈现层和 evidence 记录机制
- [ ] §2.4/§3.5 的测试全部添加并通过（§3.5 最后一条回归测试尤其重要）
- [ ] 最终报告：贴出 §4 四条命令的真实输出

## 6. 不要做的事

- 不要修改 `ConfirmAIDesignSemantics()`、`check_integration_readiness()`、`check_annotation_readiness()`、`check_rare_cell_readiness()`（如果已存在）、任何 `run_*` action 的 `prerequisites` 函数体。
- 不要改变 `RunAIAnalysis()`/`ValidateAIPlan()` 的 `status`/`valid` 字段的判定条件或取值集合。
- 不要让 `sclet_ai_format_clarification()` 具备"猜测答案"的能力——它只负责把已知的拒绝原因转成问题，不负责回答问题或建议默认值。
- 不要把 `RunBasicWorkflow()` 做成可以被 AI 调用的 action——它是给人类直接调用的确定性 facade，不进 `AIDefaultExecutionRegistry()`。
- 不要跑 `make rd`。
- 不要在仓库里留下验证脚本残留文件。

## 7. 完成后

把最终报告（§4 四条命令的真实输出）留下。审核会重点检查：`RunBasicWorkflow()` 是否真的只是顺序调用现有函数（没有重新实现任何算法逻辑）、`sclet_ai_format_clarification()` 是否真的没有新增判定能力（只做格式转换）、§3.5 最后一条回归测试是否真的证明了判定逻辑没被动过、以及 `user_decision` evidence 是否真的能被记录和查询到。
