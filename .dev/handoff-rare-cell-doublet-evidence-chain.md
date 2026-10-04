# Handoff：Rare-cell / Doublet 证据链（Phase D 第二批）

- **交接对象**：另一个 AI coding agent
- **仓库**：`/home/wang/data/source/omics/sclet`，分支 `devel`
- **前置任务状态**：Phase A/B/C 全部验收通过；Phase D 第一批（marker/DE/annotation 证据链，`check_annotation_readiness()`、`run_de_test`、`run_annotation`）已验收通过。本任务是 Phase D 第二批，跟第一批用**同一套模式**：diagnostics → action registry → evidence 登记 → 候选/确认区分。
- **不要动**：`R/ai-comparison.R`、`R/ai-evidence.R`、`R/ai-privacy.R`、`R/ai-integration-routes.R`、`R/ai-design-confirmation.R`、`run_integration`/`run_de_test`/`run_annotation` action、`check_integration_readiness()`/`check_annotation_readiness()`。这些都已验收，出问题会导致回归。

---

## 0. 一句话目标

底层的 rare-cell 检测（`RunRareCellDetection()`，density-based 方法）和 doublet 检测（`RunDoubletFinder()`，基于 `scDblFinder`）早就存在，但没有 AI-native 层。这条线的核心风险跟 annotation 那批不一样：**不是"AI 猜错了 reference"，而是"AI 只看到一个信号就把一个小群体判定为噪声删掉，结果删掉了真实但罕见的细胞类型"**。

Spec 里已经写死的规则（`.dev/ai-advanced-analysis-spec.md` 第 12.3 节 和 第 1021/1027 行）：

> 不允许仅按 cluster size 自动删除小群体；需要同时检查 QC、doublet、sample replication 和 marker；稀有群体必须报告"支持它的证据"和"替代解释"；**稀有群体必须有至少两类独立证据或明确标注低置信度**。

你要做的是把这条规则变成可执行、可验证的代码，而不是只写在文档里。

---

## 1. 硬约束（跟之前几轮一样）

```bash
# 日常验证
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'

# 完整检查，必须用 Makefile
make check
```

- `make check` 必须 `0 errors | 0 warnings | 0 notes`；
- 不要直接跑 `R CMD check`；不要跑 `make rd`（会重写 NAMESPACE，丢手工 import）；
- `R/*.R` 必须纯 ASCII；新代码注释和字符串直接写英文，不要冒非 ASCII 的风险；
- 不要泄露密钥；
- 不要破坏现有约 360 个通过的 AI 测试；
- 完成后自查并清理临时脚本：不要在仓库根目录或 `tests/testthat/` 下留下你自己验证用的草稿文件（例如 `test_audit_*.R`、`_problems/`、`*.rds` 这类）——上两轮都出现过这个问题，验证完记得删。
- 提交前跑：

```bash
cd /home/wang/data/source/omics/sclet
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"
git diff --check
```

**本机依赖现状**（实测，可以真实跑通不用 mock）：`scDblFinder`、`BiocNeighbors`、`celda`（decontX）都已安装。

---

## 2. 现有底层函数（不要重新实现，只包一层 AI-native adapter）

| 函数 | 文件 | 作用 | 已有 state |
|---|---|---|---|
| `RunRareCellDetection(sce, method="density", reduction="PCA", dims=1:20, k=20, q_threshold=0.25, rare_threshold=10, name="rareq")` | `R/priority.R:199` | 基于 mutual-kNN 密度的稀有 cluster 检测；写 `sce$rare_cluster`。**同时**调用 `sclet_set_analysis()`（供 `sclet_get_analysis(sce, "rare_cells", id)` 读取聚合结果）**和** `sclet_set_analysis_state()`（`R/priority.R:300`，标准 state record，`type="rare_cells"`，`id=name`，`method="density"`） | type `"rare_cells"`，走标准 state record 机制，跟 `doublet_finder` 是同一套（`sclet_get_state_record()`/evidence source 校验可以直接用） |
| `RunDoubletFinder(sce, ...)` | `R/data_cleaning.R:12` | 用 `scDblFinder` 计算每个细胞的 doublet score/class；写 `colData` 的 `scDblFinder.score`/`scDblFinder.class` | type `"preprocess"`，id `"doublet_finder"` |
| `RunDecontX(sce, assay="counts", ...)` | `R/data_cleaning.R:74` | ambient RNA 污染估计；写 `decontXcounts` assay 和 `colData$decontX_contamination` | **不写 analysis/state 记录**，只 `sclet_log_command()` + 加 layer——如果你要用它的输出做证据来源，需要先给它补一条 state 登记，见第 3.4 节 |

`Status()`（`R/status.R:38/84`）已经有 `has_rare_cells` health flag，读的是 `has_rare_cells(object)`（`R/analysis-accessors.R` 里定义），不需要新建这个检测逻辑，直接复用。

`sclet_state_types()`（`R/state.R:708`）已经包含 `"rare_cells"`——不需要新增 state type。`doublet_finder` 走的是标准 `preprocess` state type，也已经存在。两个都不需要动 `R/state.R`。

---

## 3. 要实现的 5 件事

### 3.1 确定性 rare-cell/doublet readiness 诊断

在 `R/ai-diagnostics.R` 里新增（跟 `check_annotation_readiness()`/`check_integration_readiness()` 同一个文件、同一种风格）：

```r
check_rare_cell_readiness(object, cluster = NULL)
```

要求：

- 检查是否已有 cluster/ident 信息（跟 `check_annotation_readiness()` 一样的检查方式）；
- 检查是否已有 PCA reduction（`RunRareCellDetection()` 默认基于 `"PCA"`）；
- 检查是否已经跑过 doublet 检测（`scDblFinder.class` 列是否存在于 `colData`）——如果没跑过，不算 blocking，但要在返回结果里标注"doublet 证据缺失，如果要做稀有群体判定，会导致独立证据数不足"；
- 返回结构参考现有两个 `check_*_readiness()` 的风格：`status`（`ready_for_diagnostic` / `not_ready`）、`checks`、`blocked_actions`、`questions`。**这里不需要 `clarification_required` 语义**——rare-cell 检测不像 integration/annotation 那样涉及"AI 猜了一个未声明的语义假设"，它更多是"输入数据是否齐备"的问题，`not_ready` 就够了。

### 3.2 新增只读聚合诊断：`summarize_small_cluster_evidence()`

这是本任务最核心的新函数，放在 `R/ai-diagnostics.R`：

```r
summarize_small_cluster_evidence(object, cluster = NULL, size_threshold = 10L)
```

功能：找出所有小于 `size_threshold` 的 cluster，对每一个聚合以下**已存在的**信号（不要重新计算，去读已经跑过的分析结果）：

- **QC 信号**：如果 `colData` 里有常见 QC 列（可以复用 `summarize_qc_by_group()` 内部的候选列名匹配逻辑，或者直接调用它，传 `group = cluster`），报告这个小群体的 QC 指标是否显著偏离全局（比如总 UMI 数、检测基因数）；
- **doublet 信号**：如果 `colData` 里有 `scDblFinder.class`，报告这个小群体里 doublet 比例；
- **marker 信号**：如果 ledger 里有已登记的 `detest`/`annotation` evidence（用 `GetAnalysisLedger()` 或 evidence registry 查），报告这个 cluster 是否有显著上调的 marker、或者是否被 `run_annotation` 标注为某个明确的细胞类型（哪怕是低置信度）；
- **sample replication 信号**：如果有 sample/batch 相关 colData 列（复用 `summarize_cluster_sample_composition()`），报告这个小群体是否分布在多个样本里，还是只出现在单一样本/批次中（只出现在一个样本里是一个警示信号，可能是批次特异的假群体而不是真实稀有细胞类型）。

返回结构：

```r
list(
    status = "available",  # 或 "not_available" + reason，跟其他诊断函数一致
    cluster = "...",
    size = <整数>,
    fraction_of_total = <数值>,
    independent_signals = list(
        qc = list(available = TRUE/FALSE, ...),
        doublet = list(available = TRUE/FALSE, ...),
        marker = list(available = TRUE/FALSE, ...),
        sample_replication = list(available = TRUE/FALSE, ...)
    ),
    n_independent_signals_available = <整数，0-4>,
    raw_values_included = FALSE
)
```

**这个函数本身不下结论**（不判定"这个群体是真实的还是噪声"），它只负责把已有的证据聚合起来，交给上层（action/evidence 那一层或者最终的人类判断）去看。这是刻意的设计：把"收集证据"和"下判断"分开，避免这个函数自己又变成一个隐藏的自动决策点。

### 3.3 Action registry：新增 `rare_cell` group

在 `R/ai-execution.R` 的 `AIDefaultExecutionRegistry()` 里，`allowed_groups` 加入 `"rare_cell"`，新增两个 action：

- `run_doublet_detection`：包 `RunDoubletFinder()`；
  - `prerequisites`：检查 `counts` assay 存在；`requireNamespace("scDblFinder")` 缺失时给可读原因；
  - `requires_confirmation = TRUE`（会写 colData）；`mutates_object = TRUE`；`allowed_state_types = "preprocess"`；`estimated_cost = "medium"`；`idempotent = FALSE`（重跑会用不同随机种子产生不同分数，不是幂等的，要如实标注，不要照抄 `run_integration` 那边的值）。

- `run_rare_cell_detection`：包 `RunRareCellDetection()`；
  - `input_schema` 至少要有 `reduction`（默认 `"PCA"`）、`rare_threshold`；
  - `prerequisites`：检查指定的 `reduction` 是否存在于 `reducedDimNames(object)`，不存在给可读原因（模仿 `run_annotation` 对缺失依赖的处理方式）；
  - `requires_confirmation = TRUE`；`mutates_object = TRUE`；`allowed_state_types = "rare_cells"`（`RunRareCellDetection()` 确实会写标准 state record，见第 2 节，跟其他 action 的写法一致，不需要特殊处理）；`estimated_cost = "medium"`；`idempotent = FALSE`（`q_threshold`/`rare_threshold` 不同会得到不同结果，且底层用了迭代式邻居标签传播，不保证跟上次完全一致）。

默认 registry（`include = "read"`）不包含这两个新 action，必须显式 `include = c("read", "rare_cell")`。

### 3.4 给 `RunDecontX()` 补一条 state 登记（如果要用它做证据来源）

如果你打算把 ambient RNA 污染估计也纳入 rare-cell 的独立证据信号之一，`RunDecontX()` 目前不写 analysis/state 记录，`sclet_ai_evidence_source()`（`R/ai-evidence.R`）要求 `source` 指向一条**状态为 `completed` 的记录**才能通过校验，所以：

- 要么给 `RunDecontX()` 也加一条 `sclet_set_analysis_state()` 登记（跟 `RunDoubletFinder()` 一样的写法），**但这会改动一个现有的公共函数，改动前确认不会破坏它现有的调用方/测试**；
- 要么本轮**不把 decontX 纳入独立证据信号**，只用 QC/doublet/marker/sample replication 这四类（这是推荐做法，风险更小，decontX 可以留到下一批）。

**建议选后者**，除非你确认改 `RunDecontX()` 是安全的。

### 3.5 Evidence 登记 + 独立证据数量门槛

执行完 `run_rare_cell_detection` 之后：

- 对每一个被判定为"小"（低于阈值）的 cluster，调用 `summarize_small_cluster_evidence()` 拿到 `independent_signals`；
- 调用 `RecordAIEvidence()` 把这个聚合结果登记为一条 evidence：
  - `kind = "deterministic_summary"`；
  - `source` 指向 `run_rare_cell_detection` 写入的 analysis id；
  - `values` 只放聚合摘要（cluster size、fraction、每类信号的可用性布尔值、以及可以安全聚合的数值比如"doublet 比例"“QC 偏离度"——不要把逐细胞值塞进去）；
  - **`claim_level` 必须根据 `n_independent_signals_available` 分级**。注意：`RecordAIEvidence()`（`R/ai-evidence.R` 第 120 行）目前只接受 `c("observed", "measured", "associated", "consistent_with")` 四个值，**不接受 `"hypothesis"`**（那是 `sclet_ai_result$findings` 那一层的 claim_level 集合，跟 evidence node 的集合不是同一个，不要搞混，调用前先确认这一点没有变化）。因此分级方式改为：
    - `n_independent_signals_available == 1`：`claim_level = "associated"`（四个允许值里最弱的一档），且 `values` 里必须显式加 `low_confidence = TRUE`；
    - `n_independent_signals_available >= 2`：`claim_level = "consistent_with"`；
    - `n_independent_signals_available == 0`：**不调用 `RecordAIEvidence()`**，直接在 action 返回值里说明"证据不足，无法判定"（比如放进 `attr(result, "sclet_ai_note")` 或者你设计的其它非 evidence 通道，只要不是登记一条 evidence 就行）。
  - **绝不允许**只凭 cluster size 本身就产生一个"这是真实稀有细胞类型"的结论性 evidence。

### 3.6 候选/确认区分：不能自动删除

跟 annotation 那批"不能覆盖 Idents()"是同一类红线：`run_rare_cell_detection` 只能**标注**（写 `sce$rare_cluster` 这个新列，这是 `RunRareCellDetection()` 已有的行为，保留），**绝不能**在这个 action 里顺手把被判定为"噪声"的细胞从对象里删除或过滤掉。删除操作如果未来要做，必须是另一个需要更高 confirmation 门槛的独立 action，本轮不做。

---

## 4. 必须补的测试

### 4.1 readiness 诊断

```r
test_that("check_rare_cell_readiness requires cluster assignment and PCA", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    result <- check_rare_cell_readiness(sce)
    expect_equal(result$status, "not_ready")
})
```

### 4.2 evidence 聚合函数的信号计数

构造一个小型 SCE，手动写入 cluster、`scDblFinder.class`、sample colData，验证 `summarize_small_cluster_evidence()` 正确统计出可用信号数量，且在信号缺失时（比如没有 doublet 列）该维度返回 `available = FALSE` 而不是报错。

### 4.3 action registry 隔离

```r
test_that("AIDefaultExecutionRegistry has rare_cell group only when explicitly included", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    default_registry <- AIDefaultExecutionRegistry(sce)
    expect_false("run_rare_cell_detection" %in% names(default_registry))
    full_registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    expect_true("run_rare_cell_detection" %in% names(full_registry))
})
```

### 4.4 端到端：真正走 `ValidateAIPlan` + `ExecuteAIPlan`

参考上一批 `run_annotation`/`run_de_test` 的端到端测试写法（`tests/testthat/test-ai-execution.R`），给 `run_rare_cell_detection` 写一条真正经过标准管线的测试，用真实小型数据（本机依赖已装，不需要 mock `BiocNeighbors`/`scDblFinder`）。

### 4.5 claim_level 分级验证——这是本任务最重要的测试

```r
test_that("rare cluster with exactly 1 independent signal is recorded as associated with low_confidence flag", {
    # 构造一个只有 cluster + 一类信号（比如只有 doublet class，没有 sample/marker）的小型 SCE
    # 跑 run_rare_cell_detection，检查产生的 evidence claim_level == "associated"
    # 且 values$low_confidence == TRUE
})

test_that("rare cluster with 2+ independent signals is recorded as consistent_with", {
    # 构造一个有 cluster + doublet class + sample 列的 SCE
    # 跑 run_rare_cell_detection，检查 claim_level == "consistent_with"
})

test_that("rare cluster with 0 independent signals does not record any evidence", {
    # 构造一个连 QC 列都没有的最小 SCE，跑 run_rare_cell_detection
    # 断言没有产生任何 ai_evidence state record（用 GetAnalysisLedger 检查）
})

test_that("rare_cluster labels never trigger automatic cell removal", {
    # 跑完 run_rare_cell_detection 后，断言 ncol(object) 不变，
    # 断言原始 counts assay 不变
})
```

---

## 5. 验收标准

- §4 五条测试类别全部有对应测试且通过；
- `summarize_small_cluster_evidence()` 是纯诊断函数，自己不下"真实/噪声"结论；
- evidence 的 `claim_level` 严格按独立信号数量分级：0 个信号不登记 evidence，1 个信号降级为 `associated` 并标记 `low_confidence = TRUE`，2 个以上才是 `consistent_with`；
- 0 个信号时不登记 evidence，只在返回值里说明证据不足；
- `run_rare_cell_detection`/`run_doublet_detection` 不删除、不过滤任何细胞，`ncol(object)` 和原始 `counts` 前后不变；
- 默认 registry 不包含这两个新 action；
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
[ FAIL 0 | WARN 0 | SKIP 1 | PASS >= 370 ]
0 errors | 0 warnings | 0 notes
```

---

## 6. 交付物

- [ ] `R/ai-diagnostics.R` 新增 `check_rare_cell_readiness()`、`summarize_small_cluster_evidence()`
- [ ] `R/ai-execution.R` 新增 `rare_cell` group、`run_doublet_detection`、`run_rare_cell_detection` 两个 action
- [ ] 对应的手写 `man/*.Rd`
- [ ] `NAMESPACE` 新增 export（不要有重复行；两个新 action 名不需要单独 export，参考 `run_integration`/`run_annotation` 的先例，只 export 诊断函数）
- [ ] `DESCRIPTION` 的 `Collate:` 如有新文件要同步（本轮如果只改现有的 `R/ai-diagnostics.R`/`R/ai-execution.R`，不需要新增条目）
- [ ] `NEWS.md` 顶部加一条本阶段条目
- [ ] `.dev/ai-advanced-analysis-spec.md` 更新：把"rare-cell / doublet diagnosis"从"还需要建设"清单里移除，在状态行和"重要状态说明"里加上做了什么、还差什么（比如"decontX 未纳入独立信号"这类范围说明）
- [ ] §4 的测试全部添加并通过
- [ ] 最终报告：贴出 §5 四条命令的真实输出
- [ ] 确认没有在仓库里留下你自己验证用的临时脚本文件

## 7. 不要做的事

- 不要只凭 cluster size 就登记一条"这是真实稀有细胞类型"的确定性 evidence。
- 不要在任何 action 里自动删除、过滤或合并细胞——本轮只做标注和证据收集，不做处理决策。
- 不要给 `summarize_small_cluster_evidence()` 加自动判断逻辑（比如内部悄悄算一个总分然后返回"是/不是稀有类型"）——它必须保持"只聚合证据、不下结论"。
- 不要动 `RunDecontX()`，除非你确认改动安全（默认按 §3.4 建议跳过它）。
- 不要动已验收的 integration/annotation/evidence/privacy/design-confirmation 模块。
- 不要跑 `make rd`。
- 不要在仓库里留下验证脚本残留文件。

## 8. 完成后

把最终报告（§5 四条命令的真实输出）留下。审核会重点检查：`summarize_small_cluster_evidence()` 是否真的没有下结论、`claim_level` 是否真的按信号数量分级（会构造一个只有 1 类信号的场景亲自验证是否降级为 `associated` 且带 `low_confidence`，以及 0 信号场景是否真的没有登记 evidence）、是否有任何路径会导致细胞被自动删除、以及是否有新测试是"绕开标准管线直接调内部函数"的假端到端测试。
