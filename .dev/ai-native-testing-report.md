# sclet R-native AI：实现与可复现测试报告

## 1. 本轮实现范围

本轮继续坚持“用户留在 R 中、AI 通过 `aisdk`、不实现外部 agent”的边界：

- 扩展 `AIDefaultExecutionRegistry()`，默认仍只开放安全的只读 action；只有显式传入 `include` 才会加入分析 action。
- 新增并登记 sclet-native action：
  - `normalize_data` → `NormalizeData()`
  - `find_variable_features` → `FindVariableFeatures()`
  - `scale_data` → `ScaleData()`
  - `run_pca` → `RunPCA()`
  - `run_umap` → `RunUMAP()`
  - `find_neighbors` → `FindNeighbors()`
  - `find_clusters` → `FindClusters()`
- action 具有输入 schema、输出 contract、前置条件、state 写入范围、成本标签和幂等性声明。
- `ValidateAIPlan()` 会投影前序 action 的声明输出，因此可以验证“先 normalize、再 HVG、再 PCA、再 graph、再 clustering”的顺序，而不会在验证阶段实际执行分析。
- 支持依赖输出绑定，例如 `${pca.output.reduction}`；引用必须在 `depends_on` 中显式声明。
- `ExecuteAIPlan()` 在每个 SCE action 返回后检查输出和 state contract；失败默认停止并登记失败记录。
- 幂等 action 支持最多 3 次 retry；`continue_on_error = TRUE` 可显式选择继续执行，并将最终状态标记为 `completed_with_errors`。
- 新增 `RunAIPlan()`，统一完成 validate → dry-run → confirmation → execute 流程。
- 修复真实 aisdk structured output 返回标量 findings/warnings 时的结果规范化问题。
- 外部 agent / Phase 4 本轮没有实现，仍按用户要求延期。

## 2. 离线、确定性测试

在仓库根目录 `/home/wang/data/source/omics/sclet` 执行：

```bash
# 聚焦本轮 AI 执行与 action catalog
Rscript -e 'devtools::test(filter = "ai-phase3", reporter = "progress")'
```

本轮结果：

```text
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 59 ]
```

覆盖内容包括：默认只读 registry、参数 schema、计划指纹、确认 token、dry-run、失败登记、顺序前置条件、native preprocessing→PCA→KNN→clustering、UMAP action、输出绑定、`RunAIPlan()`、retry 和 `continue_on_error`。

## 3. 真实 DeepSeek aisdk smoke test

`~/.Rprofile` 中的 key 只通过环境变量读取；以下命令不会输出 key：

```bash
Rscript -e 'cat("DEEPSEEK_API_KEY configured:", nzchar(Sys.getenv("DEEPSEEK_API_KEY")), "\\n")'
```

应输出 `TRUE`。真实调用示例：

```bash
Rscript - <<'RS'
pkgload::load_all('.', quiet = TRUE)
sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = matrix(c(1, 0, 3, 2, 0, 1, 4, 1, 0, 2, 1, 3), nrow = 4, ncol = 3))
)
res <- AIStatus(sce, model = "deepseek:deepseek-chat")
cat("class:", paste(class(res), collapse = ","), "\\n")
cat("task:", res$task, "\\n")
cat("structured:", isTRUE(res$metadata$structured_output), "\\n")
cat("findings:", length(res$findings), "\\n")
RS
```

本轮真实调用结果：

```text
class: sclet_ai_result,list
structured: TRUE
findings: 15
```

不要把 `DEEPSEEK_API_KEY` 写入脚本、日志、commit 或报告。若 key 未配置或网络/模型暂时不可用，离线 package tests 仍应保持可重复；在线 smoke test 属于 opt-in 检查。

## 4. 完整包检查

按仓库维护约定使用 Makefile，不要例行直接调用 `R CMD check`：

```bash
make check
```

本轮最终 `make check` 结果：

```text
Status: OK
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
Duration: 6m 47s
```

可选清理：

```bash
make clean
```

## 5. 手动复现一个受控 native workflow

下面的示例显式选择分析 action，并且仍然需要验证返回的 confirmation token：

```r
library(sclet)
library(SingleCellExperiment)

set.seed(1)
counts <- matrix(rpois(20 * 12, 5), nrow = 20, ncol = 12)
rownames(counts) <- paste0("g", seq_len(nrow(counts)))
colnames(counts) <- paste0("c", seq_len(ncol(counts)))
sce <- SingleCellExperiment(list(counts = counts))

registry <- AIDefaultExecutionRegistry(
    sce,
    include = c("read", "preprocess", "dimred", "graph", "cluster")
)
plan <- new_sclet_ai_plan(
    task = "native_chain",
    context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
    actions = list(
        list(id = "norm", action = "normalize_data"),
        list(id = "hvg", action = "find_variable_features", depends_on = "norm"),
        list(id = "scale", action = "scale_data", depends_on = "hvg"),
        list(id = "pca", action = "run_pca", params = list(ncomponents = 5), depends_on = "scale"),
        list(id = "knn", action = "find_neighbors", params = list(dims = 1:5, k = 5), depends_on = "pca"),
        list(id = "cluster", action = "find_clusters", depends_on = "knn")
    )
)

validation <- ValidateAIPlan(plan, object = sce, registry = registry)
stopifnot(isTRUE(validation$valid))

# First inspect without changing the object.
dry <- ExecuteAIPlan(sce, plan, registry, validation = validation, dry_run = TRUE)
stopifnot(identical(dry$status, "dry_run"))

# Then execute only with the validation token.
executed <- ExecuteAIPlan(
    sce, plan, registry,
    validation = validation,
    dry_run = FALSE,
    confirmation = validation$confirmation_token
)
stopifnot(identical(executed$status, "completed"))
executed$object
```

也可以使用统一入口：

```r
preview <- RunAIPlan(sce, plan, registry) # 默认 dry_run = TRUE
executed <- RunAIPlan(sce, plan, registry, dry_run = FALSE, confirm = TRUE)
```

对于包含副作用的自定义 action，必须设置 `requires_confirmation = TRUE`；不应从全局环境按字符串查找或执行函数。默认 registry 不含副作用 action。

## 6. 维护注意事项

- 修改后先跑聚焦 `devtools::test(filter = "ai-phase3")`，再跑 `make check`。
- 不要在输出中显示或保存 API key。
- `make rd` 可能重写手工维护的 `NAMESPACE`、`DESCRIPTION` roxygen 版本和无关 Rd 文件；运行后必须复核并恢复无关变化。
- `.dev/` 被 `.Rbuildignore` 排除，不会进入包构建产物。
