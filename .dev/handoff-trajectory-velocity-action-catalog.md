# Handoff：Trajectory / Velocity Action Catalog 第一批（只做 readiness + trajectory，不做 velocity/fate/spatial/multimodal 执行）

- **交接对象**：另一个 AI coding agent
- **仓库**：`/home/wang/data/source/omics/sclet`，分支 `devel`
- **前置任务状态**：Phase D 第一批（marker/DE/annotation）、第二批（rare-cell/doublet）均已验收通过（如果你在做这个任务时第二批还没做完，两者互不依赖，可以并行）。
- **不要动**：`R/ai-comparison.R`、`R/ai-evidence.R`、`R/ai-privacy.R`、`R/ai-integration-routes.R`、`R/ai-design-confirmation.R`、`run_integration`/`run_de_test`/`run_annotation`/`run_rare_cell_detection`/`run_doublet_detection` action、`check_integration_readiness()`/`check_annotation_readiness()`/`check_rare_cell_readiness()`。

---

## 0. 一句话目标，以及为什么这批范围要收窄

Spec 里把"trajectory、velocity、spatial 和 multimodal action catalog"列为一整块未完成项，但这一整块工作量和风险都远大于之前几批（annotation、rare-cell）。**这次只做其中风险最集中、最容易做错的一小块：trajectory 的 root/起点假设治理 + readiness 诊断 + 一个受限的 `run_trajectory` action。** Velocity（需要 spliced/unspliced 输入 + Python 环境）、CellRank/fate（依赖 velocity 输出）、spatial、multimodal 本轮**不做**，留给后续批次。

**核心风险**（跟之前几批是同一个模式，但换了个具体形态）：`RunSlingshot_trajectory(object, start.clus = NULL, ...)` 默认 `start.clus = NULL`——如果不指定，slingshot 会自己挑一个起点/根节点。这跟 `RunSingleR(ref=NULL)` 默认下载人类 reference、跟 AI 自称 `.design_confirmed=TRUE` 是同一类问题：**AI 不能替用户决定"轨迹从哪个 cluster 开始"这个生物学假设**，这个判断权必须留给人类。

Spec 里已经写的规则（第 10 节，`.dev/ai-advanced-analysis-spec.md`）：

> root、terminal state、reference 等关键假设必须显式记录；trajectory 不能把 cluster 顺序自动解释成时间顺序。

---

## 1. 硬约束（跟之前几批一样）

```bash
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'
make check
```

- `make check` 必须 `0 errors | 0 warnings | 0 notes`；
- 不要直接跑 `R CMD check`；不要跑 `make rd`；
- `R/*.R` 必须纯 ASCII；
- 不要破坏现有约 370+ 个通过的 AI 测试（具体数字取决于 rare-cell 那批是否已经并入，用 `devtools::test(filter="ai-")` 跑一遍看当前基线）；
- 完成后清理你自己验证用的临时脚本文件，不要留在仓库里（上几轮反复出现这个问题，这次务必注意）；
- 提交前跑：

```bash
cd /home/wang/data/source/omics/sclet
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"
git diff --check
```

**本机依赖现状**（实测）：`slingshot`、`velociraptor`、`basilisk` 都已安装。这意味着 trajectory 的端到端测试可以真实跑通、不需要 mock `slingshot`。但**本轮不涉及 velocity/RegVelo/CellRank 的实际执行**，所以不需要验证 Python 后端是否真的配置好（那是下一批的问题）。

---

## 2. 现有底层函数（本轮只包 trajectory 这一个）

| 函数 | 文件 | 作用 | 已有 state |
|---|---|---|---|
| `RunSlingshot_trajectory(object, reduction="UMAP", cluster.labels=NULL, start.clus=NULL, name="slingshot_trajectory", ...)` | `R/trajectory.R:1` | Slingshot 拟时序轨迹推断；写 `colData$slingPseudotime`；`start.clus=NULL` 时由 slingshot 自动选根 | `sclet_set_analysis()` + `sclet_set_analysis_state(type="trajectory", ...)`（`R/trajectory.R:26-40`），双重登记，跟 rare_cells 那批是同一种模式 |
| `RunSlingshot(sce, group, reduction="UMAP", start_cluster=NULL, end_cluster=NULL, reverse=FALSE, align_start=FALSE, seed=2025, name="slingshot")` | `R/runSlingshot.R:20` | 另一个 slingshot 入口，比 `RunSlingshot_trajectory` 更完整：`group` 是显式 colData 列名（不隐式读 `Idents()`）、支持 `end_cluster`、有固定默认 `seed`（确定性更明确）、`reduction` 通过 `match.arg` 限定为 `UMAP/PCA/tSNE`；跟 `RunSlingshot_trajectory` 一样走 `sclet_set_analysis()` + `sclet_set_analysis_state(type="trajectory", ...)`（`R/runSlingshot.R:81/97`）双重登记 |

`sclet_state_types()` 已包含 `"trajectory"`，不需要新增。

**本轮明确不包**（列出来是为了防止你顺手做多）：`RunVelocity()`（`R/velocity.R`，需要 velociraptor）、`RunRegVelo()`（`R/regvelo.R`，需要 basilisk Python 环境）、`RunCellRank()`/`RunCellFate()`（`R/cellrank.R`，依赖 velocity 输出）。这些留给后续批次，本轮如果顺手做了也不算数，会被要求拆出去。

---

## 3. 要实现的 4 件事

### 3.1 确定性 trajectory readiness 诊断

在 `R/ai-diagnostics.R` 新增（跟 `check_annotation_readiness()`/`check_rare_cell_readiness()` 同一种风格）：

```r
check_trajectory_readiness(object, reduction = NULL)
```

要求：

- 检查是否有可用的 reduction（`RunSlingshot_trajectory()` 默认用 `"UMAP"`，检查 `reducedDimNames(object)` 里有没有一个可用的低维嵌入，不强制要求叫 `"UMAP"`，但要报告用哪个）；
- 检查是否已有 cluster/ident 信息（trajectory 通常基于 cluster 顺序做拟时序）；
- **不检查、不建议、不推断"哪个 cluster 应该是起点"**——这是本任务最重要的边界：这个诊断函数只报告"有没有 cluster 可用"，绝不能输出类似"cluster 3 看起来像起点"这种建议。如果你想加一个"当前 cluster 分布概况"之类的只读信息帮用户自己判断，可以加，但不要包装成"推荐起点"。
- 返回结构参考现有 `check_*_readiness()` 的风格：`status`（`ready_for_diagnostic`/`not_ready`）、`checks`、`blocked_actions`。**跟 rare-cell 一样，这里不需要 `clarification_required` 语义**——真正的"起点必须显式确认"这个门槛放在下面 3.3 的 action prerequisites 里，不是在这个只读诊断里。

### 3.2 只读诊断：`summarize_trajectory_cluster_order()`（可选，视时间而定）

如果时间允许，加一个只读函数，聚合"当前 cluster 在选定 reduction 上的相对位置分布"（比如每个 cluster 在第一个主成分/嵌入维度上的均值、方差），**只报告数据分布，不做任何"这个应该是起点"的判断**。这条不是本轮强制项，如果跳过，在报告里说明。

### 3.3 Action registry：新增 `trajectory` group，`root` 必须显式确认

在 `R/ai-execution.R` 的 `AIDefaultExecutionRegistry()` 里，`allowed_groups` 加入 `"trajectory"`，新增一个 action：

```r
run_trajectory
```

**建议优先包 `RunSlingshot()`**（`R/runSlingshot.R:20`），不是 `RunSlingshot_trajectory()`：前者要求显式传 `group`（一个真实的 colData 列名，比隐式读 `Idents()` 更明确、更容易在 `prerequisites` 里校验），且有固定默认 `seed=2025`（幂等性更有保障）。如果读完代码发现两者在你的场景下有实质性差异（比如下游 evidence/state 读取方式不同、或者 `RunSlingshot_trajectory` 有 `RunSlingshot` 没有的能力），在报告里说明你的选择依据；如果时间允许两个都包也可以，但不是强制项。

**核心设计**：复用上一批 design confirmation 的思路（但不要照抄 `ConfirmAIDesignSemantics()` 那套代码，这里需求不完全一样，起点确认跟批次确认的生命周期不同）——

- `input_schema` 里 `start_cluster` 声明为**必填**参数（`required = TRUE`），类型是单个字符串或 `NULL`——但注意：如果允许 `NULL` 作为一个"合法值"传进来，那就等于没有强制要求，AI 依然可以传 `start_cluster = NULL` 混过 schema 检查。你需要在 `prerequisites` 里显式拒绝 `NULL`/缺失：

```r
prerequisites = function(object, params, planned = NULL) {
    if (is.null(params$start_cluster) || !nzchar(as.character(params$start_cluster))) {
        return("start_cluster_missing: trajectory root must be an explicit cluster label; do not rely on slingshot's automatic root selection. Inspect current cluster identities and supply start_cluster explicitly.")
    }
    idents <- Idents(object)
    if (is.null(idents)) {
        return("cluster identity is not set: run FindClusters() first")
    }
    if (!as.character(params$start_cluster) %in% unique(as.character(idents))) {
        return(paste0("start_cluster '", params$start_cluster, "' is not a current cluster identity"))
    }
    reduction <- params$reduction %||% "UMAP"
    if (!reduction %in% SingleCellExperiment::reducedDimNames(object)) {
        return(paste0("reduction '", reduction, "' not found; run the corresponding dimensionality reduction first"))
    }
    if (!requireNamespace("slingshot", quietly = TRUE)) {
        return("optional_package_missing: slingshot is required for run_trajectory; install via BiocManager::install('slingshot')")
    }
    TRUE
}
```

（上面是示意，不是要求你一字不改地照抄，但核心逻辑——`start_cluster` 缺失/NULL 时必须拒绝——不能省。）

- `requires_confirmation = TRUE`；`mutates_object = TRUE`；`allowed_state_types = "trajectory"`；`estimated_cost = "medium"`；`idempotent = TRUE`（同样输入应该得到同样输出，跟 rare-cell 那批不一样，slingshot 本身是确定性算法，不涉及随机种子——**如果你发现实际上有随机性，改成 FALSE 并说明原因**）。

- handler 执行完之后，跟之前几批一样调用 `RecordAIEvidence()` 登记：
  - `kind = "deterministic_summary"`；
  - `claim_level = "consistent_with"`（trajectory 拟时序是一种模型推断，不是直接测量）；
  - `values` 里放聚合摘要：cluster 数量、拟时序值的分布摘要（均值/分位数，不要放逐细胞的完整 pseudotime 向量）、使用的 `start_cluster`（这是一个 cluster 标签，不是逐细胞原始值，参考之前几批"group_N 匿名化"的做法，或者直接用 cluster 的原始标签——**这里要判断一下**：cluster label 本身通常不算敏感原始数据（跟 patient ID 不是一个级别），可以直接记录，但如果你觉得需要匿名化处理，参考 `run_rare_cell_detection` 的做法）；
  - **不要**把逐细胞的 pseudotime 值放进 evidence `values`（那是接近原始数据的信息，长度等于细胞数）。

默认 registry（`include = "read"`）不包含 `run_trajectory`。

### 3.4 不能把 cluster 顺序自动解释成时间顺序

跟"不能覆盖 Idents()"、"不能自动删除细胞"同类的红线：`run_trajectory` 只能写 `colData$slingPseudotime`（这是 `RunSlingshot_trajectory()` 已有行为），**不能**在 action 层面额外生成任何暗示"这是真实生物学时间"的标签或结论（比如不要自动生成"day1/day2/day3"这类时间点命名）。evidence 里如果要提到 pseudotime 分布，措辞上要清楚这是"基于选定 root 和 reduction 的相对排序"，不是绝对时间。

---

## 4. 必须补的测试

### 4.1 readiness 诊断

```r
test_that("check_trajectory_readiness requires cluster assignment and a reduction", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    result <- check_trajectory_readiness(sce)
    expect_equal(result$status, "not_ready")
})

test_that("check_trajectory_readiness never recommends a specific start cluster", {
    # 构造一个有 cluster 信息的 SCE，跑 check_trajectory_readiness，
    # 断言返回结构里不存在任何形如 "recommended_start"/"suggested_root" 的字段
})
```

### 4.2 action registry 隔离

```r
test_that("AIDefaultExecutionRegistry has trajectory group only when explicitly included", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    default_registry <- AIDefaultExecutionRegistry(sce)
    expect_false("run_trajectory" %in% names(default_registry))
    full_registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))
    expect_true("run_trajectory" %in% names(full_registry))
})
```

### 4.3 核心红线：AI 不能省略 start_cluster——真正走标准管线验证

```r
test_that("run_trajectory rejects plans without an explicit start_cluster", {
    # 构造一个有 cluster + UMAP 的小型 SCE
    # 构造一个 run_trajectory plan，params 里不传 start_cluster（或传 NULL）
    # 用 ValidateAIPlan() 验证 valid == FALSE，
    # errors 里包含 "start_cluster_missing"
})

test_that("run_trajectory succeeds end-to-end with an explicit start_cluster", {
    # 同一个对象，params 里传一个真实存在的 cluster label 作为 start_cluster
    # ValidateAIPlan() -> valid TRUE
    # ExecuteAIPlan() -> status "completed"
    # 用真实 slingshot（本机已装，不需要 mock）
})
```

### 4.4 Evidence 验证

```r
test_that("run_trajectory records consistent_with evidence without leaking per-cell pseudotime", {
    # 执行 run_trajectory 后检查 evidence claim_level == "consistent_with"
    # 且 evidence$values 里没有任何长度等于 ncol(object) 的数值向量
})
```

---

## 5. 验收标准

- §4 四条测试类别全部有对应测试且通过；
- `run_trajectory` 在缺少/为 NULL 的 `start_cluster` 时，通过标准 `ValidateAIPlan()` 管线被拒绝（不是只测 `prerequisites()` 本身，参考之前几批"假端到端测试"的教训）；
- `check_trajectory_readiness()` 不包含任何"推荐起点"的输出；
- evidence 不含逐细胞 pseudotime 向量；
- 默认 registry 不包含 `run_trajectory`；
- 确认本轮**没有**顺手实现 velocity/CellRank/fate/spatial/multimodal 的任何执行 action（只读的、纯粹为了理解现状的探索性代码读取不算，但不要往 registry 里加东西）；
- 不修改任何已验收模块；
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

## 6. 交付物

- [ ] `R/ai-diagnostics.R` 新增 `check_trajectory_readiness()`（3.2 的聚合函数可选）
- [ ] `R/ai-execution.R` 新增 `trajectory` group、`run_trajectory` action
- [ ] 对应的手写 `man/*.Rd`
- [ ] `NAMESPACE` 新增 export（诊断函数需要 export，action 名不需要，参考之前几批先例）
- [ ] `NEWS.md` 顶部加一条本阶段条目
- [ ] `.dev/ai-advanced-analysis-spec.md` 更新：把"trajectory action catalog"从"还需要建设"清单里按实际完成范围调整措辞（明确写清楚"只完成了 trajectory 的 readiness + root-confirmed 执行，velocity/CellRank/fate/spatial/multimodal 仍未开始"，不要笼统写"trajectory/velocity 已实现"）
- [ ] §4 的测试全部添加并通过
- [ ] 最终报告：贴出 §5 四条命令的真实输出，并说明 `RunSlingshot_trajectory` vs `RunSlingshot` 两者关系的调查结论（2 节表格里提到的那个需要你自己确认的点）

## 7. 不要做的事

- 不要实现 `run_velocity`/`run_cellrank`/`run_fate` 或任何 spatial/multimodal action——本轮范围明确只到 trajectory。
- 不要让 `start_cluster` 允许 `NULL` 通过 `prerequisites`。
- 不要在 `check_trajectory_readiness()` 或任何只读诊断里输出"建议的起点/根节点"。
- 不要把 pseudotime 排序包装成"真实时间点"或自动生成时间点标签。
- 不要把逐细胞 pseudotime 向量放进 evidence values。
- 不要动已验收的 integration/annotation/rare-cell/evidence/privacy/design-confirmation 模块。
- 不要跑 `make rd`。
- 不要在仓库里留下验证脚本残留文件。

## 8. 完成后

把最终报告（§5 四条命令的真实输出）留下。审核会重点检查：`start_cluster` 缺失时是否真的通过 `ValidateAIPlan()` 标准管线被拒绝（不是只测 `prerequisites()` 本身）、`check_trajectory_readiness()` 是否真的没有输出任何起点建议、evidence 是否真的不含逐细胞数据、以及是否有任何 velocity/CellRank/spatial/multimodal 的执行能力被顺手加了进来（本轮范围之外，加了要退回）。
