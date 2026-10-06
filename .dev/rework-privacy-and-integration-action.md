# 返工任务：privacy allowlist 覆盖不全 + run_integration action 端到端不可达

- **交接对象**：完成上一轮 handoff（`.dev/handoff-advanced-ai-next-phase.md`）的同一个 agent
- **仓库**：`/home/wang/data/source/omics/sclet`，分支 `devel`
- **背景**：主代理审核了 T1–T6 的实现。T2/T3/T4/T6 质量扎实，验收通过。**T1 和 T5 各有一个必须先修的问题**，本文档只覆盖这两点，不要在此基础上扩大范围。
- **范围边界**：不要重构已验收通过的部分（T2 `RunIntegrationRoutes`、T3 `CompareAIAnalyses` trade-off 检测、T4 evidence graph、T6 文档）。只改下面列出的具体点。

---

## 0. 复现方式（先跑这两段，确认问题真实存在，再动手改）

### 复现问题 1：privacy 默认值让 AI 变成瞎子

```bash
cd /home/wang/data/source/omics/sclet
Rscript -e '
pkgload::load_all(quiet=TRUE)
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, nrow=5, ncol=4)))
sce <- sclet_set_analysis_state(sce, "integration", "harmony_1", "harmony", summary=list(status="completed"))

captured <- NULL
old <- options(sclet.ai.call = function(task, context, ...) {
    captured <<- context
    list(answer="ok", findings=list())
})
on.exit(options(old))
res <- sclet:::sclet_ai_call(
    task="status",
    context = GetAnalysisLedger(sce, detail="full", include_data=TRUE),
    structured_output = FALSE
)
cat("Top-level keys actually delivered to provider:\n")
print(names(captured))
'
```

**当前输出**（错误）：

```text
[1] "schema_version" "dataset"        "capabilities"
```

**期望输出**：至少应包含 `active_view`、`analyses`、`health`、`state_records`、`workflows`、`lineage`、`warnings`、`fingerprint`（聚合状态字段本身不是隐私数据，不该被拒收）。

### 复现问题 2：`run_integration` action 通过标准执行路径永远不可达

```bash
cd /home/wang/data/source/omics/sclet
Rscript -e '
pkgload::load_all(quiet=TRUE)
set.seed(1)
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(rpois(40,5), nrow=10, ncol=4)))
SummarizedExperiment::colData(sce)$batch <- c("a","a","b","b")

registry <- AIDefaultExecutionRegistry(sce, include = c("read","integration"))
plan <- new_sclet_ai_plan(
    task = "test_plan",
    context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
    actions = list(
        list(id = "integrate", action = "run_integration",
             params = list(batch = "batch", method = "fastMNN", .design_confirmed = TRUE))
    )
)
v <- ValidateAIPlan(plan, object = sce, registry = registry)
cat("valid:", v$valid, "\n")
print(v$errors)
'
```

**当前输出**（错误）：

```text
valid: FALSE
[1] "invalid params for action integrate: unknown parameter(s): .design_confirmed"
```

**期望输出**：`valid: TRUE`（在 design 确实已被上游确认的前提下），且后续 `ExecuteAIPlan()` 能真正跑通这个 action。

---

## 1. 问题 1：修 privacy allowlist

### 1.1 根因

`R/ai-privacy.R` 第 317 行 `sclet_ai_payload_field_kind()` 的 `"aggregate"` 分类正则（第 341 行）没有覆盖以下**从 Phase 1/2 就存在的** `GetAnalysisLedger()` 顶层字段名：

```text
active_view
active_states
analyses
health
state_records
workflows
lineage
blocked_actions
quality_checks
warnings
fingerprint
dataset.modalities
```

这些字段落进 `sclet_ai_payload_field_kind()` 的默认分支 `"unknown"`，在 `standard`/`strict` 模式下被当作未分类字段拒收或要求 consent。

同时，`R/ai-adapter.R` 里 `sclet_ai_call()` 的 `enforce_privacy` 默认值是：

```r
enforce_privacy = getOption("sclet.ai.enforce_privacy", TRUE)
```

对**所有调用**（包括测试里常用的 mock）默认开启裁剪，两者叠加导致真实 payload 几乎清空。

### 1.2 要求的修法

**第一步**：扩充 `sclet_ai_payload_field_kind()`（`R/ai-privacy.R` 第 341 行那条正则），把下列**结构/状态类**字段名加入 `"aggregate"`（或视情况加入更合适的既有分类，如 `"container"`）：

```text
active_view, active_states, analyses, health, state_records,
workflows, lineage, blocked_actions, quality_checks, warnings,
fingerprint, modalities
```

**逐字段过一遍语义再决定分类**，不要简单地把整段正则改成"匹配一切"：

- 这些字段本身是**结构/状态容器**（例如 `analyses` 下面是分析记录列表，`health` 下面是布尔标志），不含原始表达矩阵、不含 cell barcode、不含患者 ID——**允许放行容器本身**，但容器内部仍然要走现有的递归 `visit()` 逐层检查（也就是说，即使 `analyses` 这个键放行了，它下面每个 record 里的具体字段仍然要过一遍 `deny`/`metadata` 检查，不能因为父键放行就整体豁免）。
- `fingerprint` 是一个已经脱敏的哈希字符串，允许放行。
- `warnings` 里如果混有自由文本（比如某条 warning message 引用了原始列名或路径），要确认现有的 `sclet_ai_payload_contains_secret_or_path()` 检查仍然在字符串层面生效（这个函数目前的调用点在 `sclet_ai_payload_atomic()` 里，只要 `warnings` 走到叶子节点仍会经过它，应该没问题，但请写测试验证）。

**第二步**：新增一条使用**真实 `GetAnalysisLedger()` 输出**（不是手写 dict）的回归测试，放在 `tests/testthat/test-ai-privacy.R`：

```r
test_that("real GetAnalysisLedger output survives standard privacy sanitization", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    sce <- sclet_set_analysis_state(
        sce, "integration", "harmony_1", "harmony",
        summary = list(status = "completed")
    )
    ledger <- GetAnalysisLedger(sce, detail = "full", include_data = TRUE)
    result <- sclet:::sclet_ai_sanitize_context(NULL, ledger, privacy = "standard")

    must_survive <- c(
        "schema_version", "dataset", "active_view", "analyses",
        "health", "state_records", "workflows", "lineage",
        "warnings", "fingerprint", "capabilities"
    )
    missing <- setdiff(must_survive, names(result$payload))
    expect_length(missing, 0L)

    # matrix values and raw identifiers must still be absent
    serialized <- paste(capture.output(str(result$payload)), collapse = " ")
    expect_false(grepl("^[0-9]+x[0-9]+ matrix", serialized))
})
```

**第三步**：再补一条端到端测试，验证 `sclet_ai_call()` 在**默认设置**下（不显式传 `enforce_privacy`）用真实 ledger 调用时，模型收到的 payload 里保留了必要字段：

```r
test_that("sclet_ai_call default privacy still delivers ledger structure to the provider", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    sce <- sclet_set_analysis_state(
        sce, "integration", "harmony_1", "harmony",
        summary = list(status = "completed")
    )
    captured <- NULL
    old <- options(sclet.ai.call = function(task, context, ...) {
        captured <<- context
        list(answer = "ok", findings = list())
    })
    on.exit(options(old), add = TRUE)

    sclet:::sclet_ai_call(
        task = "status",
        context = GetAnalysisLedger(sce, detail = "full", include_data = TRUE),
        structured_output = FALSE
    )
    expect_true(all(c("analyses", "health", "active_view") %in% names(captured)))
})
```

### 1.3 验收标准

- 上面两条新测试通过；
- §0 "复现问题 1" 的命令重跑后，输出里包含 `analyses`、`health`、`active_view`、`state_records`、`workflows`、`lineage`、`warnings`、`fingerprint`；
- 仍然要保证原有隐私测试全部通过（尤其是 `test-ai-privacy.R` 里已有的 secret/matrix/raw-identifier 拒收测试）——**不能为了放行结构字段而放宽对原始值/密钥/路径的拒收**；
- `test-ai-integration.R` 第 125/148 行那两个 `AIStatus()` 测试重跑后，**不应再出现 "outbound payload requires user consent" 和 "Redacted xxx: field is not allowlisted" 的 WARNING**（除非确实存在需要脱敏的字段，那种情况请在报告里说明是哪个字段、为什么该脱敏）。

---

## 2. 问题 2：修 `run_integration` 无法通过标准执行路径

### 2.1 根因

`R/ai-execution.R` 第 472 行 `run_integration` action：

- `prerequisites`（第 500 行起）要求 `isTRUE(params$.design_confirmed)`；
- 但 `input_schema`（第 491 行）里**没有声明** `.design_confirmed`；
- `ValidateAIPlan()` 的参数校验逻辑（`sclet_ai_validate_action_params`，约在本文件第 96–159 行）会把任何不在 `input_schema` 里的键判定为 `"unknown parameter(s)"` 并拒绝整个 plan。

结果是：这个 action 的前置条件**只能通过一个必然被参数校验拒绝的参数来满足**——形成了一个自己无法通过的死循环。当前测试（`tests/testthat/test-ai-execution.R` 第 24/37/50/64 行）只是直接调用 `action$prerequisites(sce, list(...))`，绕开了 `ValidateAIPlan()`，所以测试是绿的，但真实调用链（`AIPlanAnalysis` → `ValidateAIPlan` → `ExecuteAIPlan`）打不通。

### 2.2 要求的修法（任选一种，但要在报告里说明为什么选它）

**方案 A（更简单，但要接受局限）**：把 `.design_confirmed` 加入 `run_integration` 的 `input_schema`，类型声明为 `logical`：

```r
input_schema = list(
    batch = list(type = "character", required = TRUE),
    method = c("fastMNN", "Harmony", "scVI"),
    name = "character",
    dims = "integer_vector",
    features = "character_vector",
    layer = "character",
    reduction = "character",
    .design_confirmed = "logical"
),
```

**局限**（必须在报告里承认）：这样一改，任何能生成 plan 的调用方（包括 AI 自己生成的 plan JSON）都可以自称 `.design_confirmed = TRUE`，`prerequisites` 无法验证这个声明是否真的经过了人工确认。这是"能跑通"和"真的安全"之间的权衡，不是白板修复。

**方案 B（更稳，工作量更大）**：不接受 plan 里任意声称的布尔标志，而是让 `prerequisites` 去查**已经写入 ledger 的确认记录**——例如复用 `check_integration_readiness()` 或 `RunIntegrationRoutes()` 里已经落地的「design 已确认」状态（比如查 `sclet_get_state_record(object, "integration", <route_id>)$summary$design_semantics == "confirmed"`，或者定义一个新的、专门记录"用户已确认这批 design 语义"的 state 记录，`prerequisites` 去查这个记录是否存在且与当前 object fingerprint 匹配）。这种方式下，`.design_confirmed` **不需要**出现在 `input_schema` 里，因为它不再是一个可由调用方随意声明的参数，而是从对象状态里读出来的事实。

**如果时间有限，先做方案 A 让链路跑通，但必须在报告里写清楚方案 A 的局限，并把方案 B 列为后续 TODO 写进 `.dev/ai-advanced-analysis-spec.md`**（不要不写就假装问题已经彻底解决）。

### 2.3 必须补的测试：真正走 `ValidateAIPlan` + `ExecuteAIPlan` 的端到端用例

在 `tests/testthat/test-ai-execution.R`（或新建 `test-ai-execution-integration.R`）里加一条**不绕开标准管线**的测试，参考骨架：

```r
test_that("run_integration action can be validated and executed through the standard plan pipeline", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    rownames(sce) <- paste0("g", seq_len(10))
    colnames(sce) <- paste0("c", seq_len(4))
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 5)
    sce <- RunPCA(sce, ncomponents = 2)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- new_sclet_ai_plan(
        task = "test_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(
                id = "integrate", action = "run_integration",
                params = list(batch = "batch", method = "fastMNN", .design_confirmed = TRUE)
            )
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)

    # mock the real integration call so the test does not depend on batchelor internals
    executed <- testthat::with_mocked_bindings({
        ExecuteAIPlan(
            sce, plan, registry,
            validation = validation, dry_run = FALSE,
            confirmation = validation$confirmation_token
        )
    }, .package = "sclet", RunIntegration = function(object, ...) {
        args <- list(...)
        SingleCellExperiment::reducedDim(object, "correctedPCA") <- matrix(1:8, ncol = 2)
        sclet_set_analysis_state(
            object, "integration", args$name %||% "fastmnn",
            method = "mocked_integration", inputs = list(batch = args$batch),
            active = FALSE
        )
    })
    expect_equal(executed$status, "completed")
})
```

如果你选了方案 B，这条测试的 `params` 里就不会有 `.design_confirmed`，而是要先跑一次能落地"design 已确认"状态的步骤（比如先跑 `RunIntegrationRoutes()` 的部分逻辑，或直接写入你新定义的确认记录），再验证 `run_integration` 的 `prerequisites` 能读到它——请按你实际的实现调整这条测试，但**核心要求不变**：必须经过 `ValidateAIPlan()` 判定 `valid = TRUE`，并且 `ExecuteAIPlan()` 能真正跑到 `status = "completed"`。

### 2.4 验收标准

- §0 "复现问题 2" 的命令重跑后，`valid: TRUE`；
- 新增的端到端测试（2.3）通过，且**没有绕开 `ValidateAIPlan`/`ExecuteAIPlan`**；
- 如果选方案 A：报告里写清楚局限，并在 Spec 里记一条方案 B 的 TODO；
- 如果选方案 B：`.design_confirmed` 不再需要出现在 `input_schema`，且有测试证明"未确认" 状态下 `ValidateAIPlan` 仍然拒绝该 action。

---

## 3. 提交前必须跑的命令（照抄，不要自己删减）

```bash
cd /home/wang/data/source/omics/sclet

# 1. 纯 ASCII
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"

# 2. 无空白错误
git diff --check

# 3. 复现命令必须显示问题已修复（贴真实输出，不要只写"已修复"）
Rscript -e '... 复现问题 1 的脚本 ...'
Rscript -e '... 复现问题 2 的脚本 ...'

# 4. 聚焦 AI 测试全绿，且不应再出现本文档 §1.3 提到的那两条 WARNING
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'

# 5. 完整包检查干净
make check
```

## 4. 交付物

- [x] `R/ai-privacy.R` allowlist 修复
- [x] `R/ai-execution.R`（和/或新文件，如选方案 B）里 `run_integration` 的修复
- [x] 新增的三条测试（§1.2 两条 + §2.3 一条），全部通过
- [x] 若选方案 A：`.dev/ai-advanced-analysis-spec.md` 里补一条方案 B 的 TODO
- [x] 最终报告：贴出 §3 五条命令的真实输出，并说明选了哪个方案、为什么

## Status synchronization (2026-10-06)

The allowlist and integration confirmation path are implemented in the current `devel` checkout. The completed path uses design-confirmation scheme B; the scheme-A TODO item is not applicable. The shared verification baseline is `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 1000 ]` for the focused AI tests and `0 errors | 0 warnings | 0 notes` for `make check`.

## 5. 不要做的事

- 不要重新设计 `RunIntegrationRoutes()`、`CompareAIAnalyses()`、evidence graph——这些已经验收通过。
- 不要为了让测试通过而放宽 `sclet_ai_evidence_bad_field()` 或 payload 里对 secret/matrix/raw-identifier 的拒收规则。
- 不要把 `enforce_privacy` 默认值改回 `FALSE`——默认严格是对的，问题是 allowlist 覆盖不全，不是"要不要默认开启"这件事本身。
- 不要只改 `R/ai-privacy.R` 的正则就假装测试过了——必须用真实 `GetAnalysisLedger()` 输出验证，手写 dict 不算数（这正是上一轮出问题的原因）。
