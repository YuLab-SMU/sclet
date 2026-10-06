# Handoff：关闭 `.design_confirmed` 自签漏洞（方案 B）

- **交接对象**：完成前两轮任务的同一个 agent
- **仓库**：`/home/wang/data/source/omics/sclet`，分支 `devel`
- **背景**：上一轮返工（`.dev/rework-privacy-and-integration-action.md`）用方案 A 让 `run_integration` action 通过了标准执行管线，但方案 A 的局限是 `.design_confirmed` 是 plan 里任意一个可以被自称为 `TRUE` 的布尔参数——**AI 生成的 plan 自己就能签发这张"已确认"通行证**，`prerequisites` 完全无法验证这个声明是否真实。这次要把这个洞堵上。
- **已验收、不要动的部分**：`RunIntegrationRoutes()`、`CompareAIAnalyses()`、evidence graph（`RecordAIEvidence`/`ValidateAIEvidenceRefs`/`sclet_ai_evidence_independence`）、privacy allowlist、T6 文档同步。只改 `run_integration` action 的 design 确认机制。

---

## 0. 复现当前漏洞（先跑，确认问题真实存在）

```bash
cd /home/wang/data/source/omics/sclet
Rscript -e '
pkgload::load_all(quiet=TRUE)
set.seed(1)
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(rpois(40,5), nrow=10, ncol=4)))
SummarizedExperiment::colData(sce)$batch <- c("a","a","b","b")

registry <- AIDefaultExecutionRegistry(sce, include = c("read","integration"))
# 攻击场景：AI 自己生成的 plan，凭空声称设计已确认，从未经过任何用户交互
plan <- new_sclet_ai_plan(
    task = "ai_self_generated_plan",
    context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
    actions = list(
        list(id = "integrate", action = "run_integration",
             params = list(batch = "batch", method = "fastMNN", .design_confirmed = TRUE))
    )
)
v <- ValidateAIPlan(plan, object = sce, registry = registry)
cat("valid (应该是 FALSE，因为没有任何真实确认记录):", v$valid, "\n")
'
```

**当前输出**（问题）：`valid: TRUE` —— 没有人确认过 `batch` 是真实技术批次，plan 照样通过验证。这就是要修的漏洞。

---

## 1. 要求的修法

### 1.1 设计

不再让 `run_integration` 的 `prerequisites` 相信 plan 参数里的自称声明，而是让它去查**已经写入 ledger 的、独立于当前 plan 的确认事实**。

`R/ai-integration-routes.R` 里 `RunIntegrationRoutes()` 已经在每条路线的 state record 里写了 `summary$design_semantics == "confirmed"`（见第 118–178 行 `sclet_ai_register_baseline_route()` / `sclet_ai_run_one_route()`）。但这只覆盖了"已经跑过一次 route"的场景，`run_integration` 这个更底层的单步 action 需要一个**跑之前就能查的确认记录**，不能依赖"必须先跑过一次才算确认"这种循环依赖。

因此新增一个专门的确认记录类型，而不是复用 `integration` state type：

**第一步**：在 `R/state.R` 的 `sclet_state_types()`（约第 686–711 行）里新增一个类型：

```r
"ai_design_confirmation"
```

（参考上一轮加 `"ai_evidence"` 的方式，直接在向量里追加一项，不要动其他项。）

**第二步**：新增一个用户/上游显式调用的确认函数，建议命名 `ConfirmAIDesignSemantics()`（放在 `R/ai-execution.R` 或新建 `R/ai-design-confirmation.R`，你可以自行决定文件归属，但要同步进 `DESCRIPTION` 的 `Collate:`）：

```r
ConfirmAIDesignSemantics <- function(object, design) {
    # design 至少要有 $batch；可选 $condition/$subject/$sample
    # 校验 design 里声明的每个列名都真实存在于 colData(object)
    # 把确认事实写入 state type "ai_design_confirmation"，
    #   id 建议用 design 的规范化摘要（比如 batch 列名 + condition 列名拼出的稳定 id），
    #   summary 至少包含: design（哪些角色对应哪些列名）、
    #                     object_fingerprint（GetAnalysisLedger(object)$fingerprint）、
    #                     confirmed_at
    # 返回更新后的 object
}
```

这个函数**不属于 AI action registry**，不应该被 AI 自动调用——它代表的是"人类/上游调用方明确知道 batch/condition 语义并主动确认"这个事实，必须由外部显式调用触发，AI 计划本身不能生成对它的调用。

**第三步**：改 `run_integration` 的 `prerequisites`（`R/ai-execution.R` 第 500 行左右）：

- 从 `input_schema` 里**删掉** `.design_confirmed`（这是方案 A 留下的洞，必须删）；
- `prerequisites` 改成去查 `sclet_get_state_records(object, "ai_design_confirmation")`，找一条：
  - `summary$design$batch` 等于本次 `params$batch`；
  - `summary$object_fingerprint` 等于当前 `GetAnalysisLedger(object)$fingerprint`（**这一点很重要**：如果对象在确认之后又发生了变化——比如换了别的 batch 列语义——旧的确认记录不应该继续有效）；
  - 如果找不到匹配记录，返回明确原因，例如：

    ```r
    return(paste0(
        "design_semantics_not_confirmed: call ConfirmAIDesignSemantics(object, design = list(batch = '",
        batch, "', ...)) before running integration"
    ))
    ```

### 1.2 为什么要绑定 fingerprint

如果不绑定 fingerprint，会出现这种攻击面：用户对数据集 A 确认了 `batch` 语义，然后对象被别的操作换了列内容（比如同一个 `batch` 列名下的值语义变了），旧确认记录依然"有效"，AI 就能拿旧的许可证在新语义下瞎跑。绑定 fingerprint 之后，任何会改变 ledger fingerprint 的操作都会让旧确认失效，必须重新确认。

### 1.3 与现有 `RunIntegrationRoutes()` 的关系

`RunIntegrationRoutes()` 的调用契约里，`design` 是**调用方直接传入**的（不经过 AI plan），所以它本来就不受这个漏洞影响，**不需要改**。但为了保持概念一致，建议（非强制）让 `RunIntegrationRoutes()` 在成功执行后，也顺手调用一次 `ConfirmAIDesignSemantics()` 把确认记录写下来，这样如果用户后续想通过 `run_integration` action 对同一批 design 语义做更细粒度的单步操作，就不需要重复确认。这一点如果时间不够可以跳过，不是本次强制验收项。

---

## 2. 必须补的测试

### 2.1 负面测试：自称确认不再生效

```r
test_that("run_integration rejects plans that only self-declare design confirmation", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "attacker_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")
            # 注意：这里不再传 .design_confirmed，因为它已经从 input_schema 里删除
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("design_semantics_not_confirmed", v$errors)))
})
```

### 2.2 正面测试：真实确认后才能通过

```r
test_that("run_integration succeeds only after ConfirmAIDesignSemantics is called", {
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

    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "confirmed_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid)

    executed <- testthat::with_mocked_bindings({
        ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE,
            confirmation = v$confirmation_token)
    }, .package = "sclet", RunIntegration = function(object, ...) {
        args <- list(...)
        SingleCellExperiment::reducedDim(object, "correctedPCA") <- matrix(1:8, ncol = 2)
        sclet:::sclet_set_analysis_state(object, "integration", args$name %||% "fastmnn",
            method = "mocked_integration", inputs = list(batch = args$batch), active = FALSE)
    })
    expect_equal(executed$status, "completed")
})
```

### 2.3 负面测试：确认记录随 fingerprint 失效

```r
test_that("stale design confirmation (fingerprint mismatch) is rejected", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))

    # 对象发生变化，fingerprint 改变（比如新增一个分析记录）
    sce <- sclet:::sclet_set_analysis_state(
        sce, "preprocess", "norm_1", "NormalizeData", summary = list(status = "completed")
    )

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "stale_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
})
```

> 如果第三条测试因为 fingerprint 计算方式的细节导致不好复现（比如某些操作不改变 fingerprint），可以换一种方式制造"确认记录不匹配当前对象"的场景，但**必须保留"确认记录会过期"这个核心断言**，不要删掉这条验收点。

---

## 3. 验收标准

- §0 的复现脚本重跑后，`valid: FALSE`；
- §2 三条测试全部通过；
- 旧的 `test-ai-execution.R` 里"通过标准管线执行 `run_integration`"的端到端测试（上一轮加的）要相应改写，改成先调用 `ConfirmAIDesignSemantics()`，不再在 plan 参数里塞 `.design_confirmed`；
- `input_schema` 里确认已删除 `.design_confirmed`；
- `sclet_state_types()` 新增 `"ai_design_confirmation"`；
- `.dev/ai-advanced-analysis-spec.md` 里把"方案 B TODO"那一行更新为"已实现"，并简要说明新增的 `ConfirmAIDesignSemantics()` API；
- 完整跑一遍：

```bash
cd /home/wang/data/source/omics/sclet
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"
git diff --check
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'
make check
```

期望：

```text
[ FAIL 0 | WARN 0 | SKIP 1 | PASS >= 331 ]
0 errors | 0 warnings | 0 notes
```

---

## Status synchronization (2026-10-06)

The design-confirmation hardening is implemented in the current `devel` checkout: `ConfirmAIDesignSemantics()`, the `ai_design_confirmation` state, stale fingerprint/value checks, removal of the AI-self-declared `.design_confirmed` input, and the corresponding tests are present. The shared verification baseline is `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 1000 ]` for the focused AI tests and `0 errors | 0 warnings | 0 notes` for `make check`. Any `design_confirmed = TRUE` text remaining in route summaries is validation metadata, not an AI-controllable plan parameter.

## 4. 不要做的事

- 不要重新设计 `RunIntegrationRoutes()`、`CompareAIAnalyses()`、evidence graph、privacy allowlist——这些已经验收通过。
- 不要把 `ConfirmAIDesignSemantics()` 注册进 `AIDefaultExecutionRegistry()` 的任何 action 里——它是给人类/上游显式调用的边界函数，不是 AI 可以自己触发的 action。如果你觉得需要一个只读的"检查是否已确认"的诊断 action，可以加，但**确认动作本身不能由 AI 触发**。
- 不要只做 fingerprint 之外的宽松匹配（比如只匹配 `batch` 列名，不检查 fingerprint）——那样等于没有解决"对象变了但确认没变"的问题。
- 不要为了让旧测试通过而保留 `.design_confirmed` 参数——必须真的从 `input_schema` 删除，否则又是一个可以被绕过的后门。

## 5. 完成后

把最终报告（§3 四条命令的真实输出）留下。主代理会重新用 §0 的攻击脚本验证一次，并检查 `ConfirmAIDesignSemantics()` 是否确实不在任何 action registry 里可达。
