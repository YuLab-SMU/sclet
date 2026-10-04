# Handoff：sclet 高级 AI 分析 —— 下一阶段实现任务

- **交接对象**：另一个 AI coding agent
- **仓库**：`/home/wang/data/source/omics/sclet`
- **分支**：`devel`（直接在此分支工作，不要切分支）
- **任务性质**：在已完成的只读骨架上，实现**真实 integration 路线执行 + evidence graph 整合 + strict privacy 收口**
- **完成后**：由原 agent 审核，你只需交付「可复现的代码 + 测试 + 文档 + 验收证据」

---

## 0. 一句话目标

当前 sclet AI 已经能**画像、诊断、登记 evidence、裁剪出站 payload、比较已记录的 integration 摘要**，但**还不能真正执行一条高级分析路线**。

你要做的是：让 AI 在**用户确认实验设计**之后，能够**真实执行** integration（fastMNN / Harmony / scVI）路线，并把结果登记为**可比较、可验证、带证据链**的 analysis 记录，从而让 `CompareAIAnalyses()` 从「比较已有摘要」升级为「比较 AI 真实产生的多路线结果」。

---

## 1. 硬约束（违反会被打回）

### 1.1 构建与检查

```bash
# 日常验证：聚焦测试
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'

# 完整检查：必须用 Makefile，不要直接 R CMD check
make check
```

- **必须** `make check` 通过：`0 errors | 0 warnings | 0 notes`
- **不要**直接调用 `R CMD check` / `rcmdcheck`
- 运行 `make rd` 会**重写 NAMESPACE 并可能丢掉手工维护的 import**（曾导致 30+ 无关测试失败）。**默认不要跑 `make rd`**；改为手工维护 `NAMESPACE` + 手写 `man/*.Rd` + 手工在 `DESCRIPTION` 的 `Collate:` 里加文件名。

### 1.2 R 源码必须纯 ASCII

`R CMD check` 会因非 ASCII 报 WARNING。中文注释/字符串一律用 `\uXXXX` escape，例如：

```r
# 正确
grepl("\\u6279\\u6b21", question, ignore.case = TRUE)
```

**不要**在 `R/*.R` 里直接写中文。写完自查：

```bash
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R
```

### 1.3 不要在输出中泄露密钥

**绝对不要**打印、写入、提交 `DEEPSEEK_API_KEY` 或任何 token。

### 1.4 在线测试双重门控（必须保持）

只在下面**同时成立**时才访问网络：

```r
nzchar(Sys.getenv("DEEPSEEK_API_KEY")) &&
identical(Sys.getenv("SCLET_RUN_ONLINE_TESTS"), "true")
```

`make check` 时在线测试必须 skip，不允许联网。

### 1.5 兼容性

不得破坏现有 223 个通过的 AI 测试。现有公开 API 签名不要改（可以**新增可选参数**）。

---

## 2. 当前已实现的状态（不要重复实现）

### 2.1 模块与文件

| 文件 | 作用 |
|---|---|
| `R/ai-context.R` | `GetAnalysisLedger()` bounded ledger |
| `R/ai-adapter.R` | `sclet_ai_call()` provider adapter + privacy 接入 |
| `R/ai-result.R` | `sclet_ai_result` 契约与校验 |
| `R/ai-tools.R` | 只读 AI 工具 |
| `R/ai-functions.R` | `AskAI()`、`RunAIAnalysis()` |
| `R/ai-planning.R` | `AIPlanAnalysis()`、`ValidateAIPlan()` |
| `R/ai-execution.R` | `AIAction()`、`AIDefaultExecutionRegistry()`、`ExecuteAIPlan()` |
| `R/ai-profile.R` | `GetAIProfile()` 数据集画像 |
| `R/ai-diagnostics.R` | 确定性诊断（QC/组成/PCA 关联/小 cluster/readiness） |
| `R/ai-investigate.R` | `AIInvestigate()` 只读调查 |
| `R/ai-evidence.R` | `RecordAIEvidence()`、`ValidateAIEvidenceRefs()` |
| `R/ai-privacy.R` | `sclet_ai_payload_policy()`、`sclet_ai_sanitize_context()` |
| `R/ai-comparison.R` | `CompareAIAnalyses()` 只读路线比较 |

测试：`tests/testthat/test-ai-{beginner,comparison,diagnostics,evidence,integration,investigate,online,phase3,privacy,profile}.R`
文档：`man/*.Rd`（手写），`.dev/ai-advanced-analysis-spec.md`（Spec，中文）

### 2.2 关键 API 契约

**数据集画像：**

```r
GetAIProfile(object, max_sparse_elements = 1000000L)
# 返回 schema_version/dataset/design/qc/structure/analysis_state/
#      diagnostics/capabilities/cost_estimates/privacy/fingerprint
# 绝不返回原始矩阵、barcode、metadata 原始值
```

**确定性诊断：**

```r
summarize_qc_by_group(object, group)
summarize_pca_metadata_association(object, metadata, reduction = NULL)
summarize_cluster_sample_composition(object, sample, cluster = NULL)
summarize_small_clusters(object, threshold = 10L, cluster = NULL)
check_integration_readiness(object, design = NULL)
# 缺失输入 -> list(status = "not_available", reason = ...)
```

**只读调查：**

```r
AIInvestigate(object, question, design = NULL,
              scope = c("diagnosis", "advanced", "report"),
              model = NULL, structured_output = TRUE, ...)
# 不执行 action；缺 design 语义时返回 clarification_required
```

**Evidence：**

```r
RecordAIEvidence(object, evidence, source = NULL, parents = NULL, scope = NULL)
ValidateAIEvidenceRefs(object, refs, fingerprint = NULL)
# evidence$kind ∈ deterministic_summary|plot|test|state|user_decision
# claim_level ∈ observed|measured|associated|consistent_with
# 值必须是聚合值；字符串只允许 ^(group|cluster|condition|sample|route)_[0-9]+$
# 写入 state type "ai_evidence"（已在 sclet_state_types() 注册）
```

**Privacy：**

```r
sclet_ai_payload_policy(object, context, privacy = c("strict","standard","local"))
sclet_ai_sanitize_context(object, context, privacy = ..., user_consent = FALSE)
# 返回 allowed_fields/redacted_fields/aggregation_threshold(10L)/
#      estimated_tokens/requires_user_consent/warnings/payload_fingerprint
# sclet_ai_call() 已接入；mock 测试路径默认跳过，
# 可用 options(sclet.ai.enforce_privacy = TRUE) 强制启用
```

**路线比较（当前是只读骨架）：**

```r
CompareAIAnalyses(object, ids = NULL,
                  criterion = c("evidence","stability","biological_preservation","cost"))
# 无记录 -> status = "not_available", reason = "no_comparable_integration_records"
# 有记录 -> status = "available", routes/baseline/metrics/tradeoffs/
#           recommendation = NULL, execution = list(allowed=FALSE, performed=FALSE)
# 只读 ledger 已有摘要，不执行 integration
```

**Action registry：**

```r
AIAction(name, handler, description = NULL,
         prerequisites = function(object, params) TRUE,
         returns = c("sce","value"), requires_confirmation = TRUE,
         input_schema = list(), output_schema = list(),
         mutates_object = NULL, allowed_state_types = character(),
         estimated_cost = "low", idempotent = FALSE)

AIDefaultExecutionRegistry(object, include = "read")
# include ∈ read | preprocess | dimred | graph | cluster | all
# 默认只读；分析 action 必须显式 opt-in
```

### 2.3 依赖现状（**重要**）

在本机实测：

```text
harmony    : NOT installed   ← 不能假设可用
batchelor  : installed       ← fastMNN 可用
basilisk   : installed       ← scVI 走 basilisk
scvi-tools : basilisk 环境内，未验证
```

`RunIntegration(object, method = c("fastMNN","Harmony","scVI"), batch, features, layer, reduction, name, ...)` 已存在（`R/integration.R`），会写 analysis state（type `integration`）。

**因此：所有测试必须 mock，不能依赖 harmony 真装。**

---

## 3. 你要实现的 6 个任务

> 建议按顺序做，每个任务独立可验证。不要一次性大改。

### T1. Integration action 注册（走 registry，不绕过确认）

**目标**：让 AI 能在 plan 里提出 integration 路线，并经过既有 validation / dry-run / confirmation 边界。

**做法**：

- 在 `AIDefaultExecutionRegistry()` 增加新 group，例如 `"integration"`（并加入 `allowed_groups`，同时让 `all` 包含它）。
- 新 action 建议：
  - `run_integration_fastmnn`
  - `run_integration_harmony`
  - `run_integration_scvi`
  - 或统一 `run_integration` + `method` 枚举参数（二选一，但要在文档里说清）
- 每个 action 必须：
  - `input_schema` 声明 `batch`（必填）、`method`、`name`、`dims` 等；
  - `prerequisites` 检查：
    - object 有可用的 `batch` 列；
    - **`design` 未确认时返回失败原因**（不得靠猜 batch）；
    - 依赖包缺失时返回可读原因（harmony 未装 → `FALSE` + 原因，**不是**抛原始错误）；
  - `allowed_state_types = "integration"`；
  - `mutates_object = TRUE`；
  - `requires_confirmation = TRUE`（涉及修改对象）；
  - `estimated_cost = "medium"` 或 `"high"`（scVI）；
  - `idempotent = FALSE`（除非你论证清楚）。
- **不得**把 harmony 变成硬依赖（保持 Suggests + `requireNamespace` 守卫）。

**验收**：

- 新增测试证明：未确认 design → plan 校验失败或 `clarification_required`；
- harmony 未安装时 action 返回**可读的 typed 原因**，不抛裸错误；
- 默认 registry（`include = "read"`）**不包含** integration action。

---

### T2. 真实多路线执行与结果登记

**目标**：能真实跑出 ≥2 条路线（例如 raw baseline + fastMNN），并把结果登记成可比较记录。

**做法**：

- 新增一个受控入口（建议名 `RunIntegrationRoutes()`，具体命名你可定，但要写进文档和 Spec）：
  - 输入：object、`design`（必须显式，含 batch/condition）、`routes`（例如 `c("raw","fastMNN","Harmony")`）、`confirm`；
  - 每条路线执行后写入 analysis state（type `integration`），并在 `summary` 中写入**可比较指标**；
  - 返回结构和 `RunAIAnalysis()` 保持一致的风格（`object`/`plan`/`validation`/`preview`/`execution`/`status`/`report` 之类的自洽结构）；
  - `confirm = "ask"` 只在交互式提示；非交互式返回 dry-run 不执行（与 `RunAIAnalysis()` 行为一致）。
- **保护原始对象**：
  - 不得覆盖原始 counts；
  - 每条路线写**独立命名**的 reduction（例如 `fastMNN` / `Harmony` / `scVI`）；
  - 保留 raw baseline。
- **指标契约**（每条路线必须登记，缺就写 `not_available` + 原因，**不许猜**）：

```r
list(
    name = "graph_lisi",              # 或其他明确定义的指标
    value = 0.82,
    direction = "higher_is_better",   # higher_is_better|lower_is_better|target_range
    scope = list(samples = ..., groups = ...),
    method = "...", parameters = list(...),
    baseline_ref = "raw",
    uncertainty = list(status = "not_reported"),
    evidence_id = "evidence_..."
)
```

至少覆盖：

- batch mixing（如 LISI / kBET，或你明确选定并写清方法的指标）
- biological preservation（condition/label separation、marker/program 保留）
- cluster stability（重采样或参数扰动下的 ARI/NMI 等）
- runtime / 成本

**验收**：

- 至少 2 条路线的端到端测试（**mock 掉真实 integration 调用**，用 `testthat::with_mocked_bindings(..., .package = "sclet")` 或 mock `RunIntegration`）；
- 测试证明原始 assay / reduction 未被覆盖；
- 测试证明每条路线都有指标或显式 `not_available`；
- 测试证明**不允许把单一分数当作最佳路线**。

---

### T3. `CompareAIAnalyses()` 升级为真实路线比较

**目标**：让比较函数读懂 T2 产生的记录，输出可判读的 trade-off，而不是只回显摘要。

**做法**：

- 保持现有返回字段兼容（`status`/`criterion`/`routes`/`baseline`/`metrics`/`tradeoffs`/`recommendation`/`execution`/`fingerprints`/`warnings`）。
- 新增：
  - **按 criterion 排序/分组**（`evidence`/`stability`/`biological_preservation`/`cost`）；
  - **trade-off 检测**：例如 batch mixing 改善但 biological preservation 下降 → 必须显式标记为 trade-off，**不得**自动选优；
  - `recommendation` 保持 `NULL`，或只有在用户显式请求时才给出**带不确定性说明**的候选（默认仍为 `NULL` 最安全）；
  - 指标缺失 → `not_available`，不猜。
- 仍然**只读**：`execution = list(allowed = FALSE, performed = FALSE)`。

**验收**：

- 两个方向冲突的模拟路线 → 输出 trade-off 标记，且 `recommendation` 仍为 `NULL`；
- 指标缺失 → `not_available` 而非编造；
- 现有 `test-ai-comparison.R` 全部仍通过。

---

### T4. Evidence graph 与 lineage 连接

**目标**：把 evidence node 真正接进 analysis lineage，而不是孤立记录。

**做法**：

- evidence node 增加/填充 `parents`、`dependency_group`；
- 提供 lineage 查询，能回答：这条 claim 由哪些 analysis / evidence 支撑，是否共享同一原始输入；
- **证据独立性判定由 R 计算**（共享 raw input / 同一聚类分组 / 同一 reference 不算独立），AI 只解释结果；
- 校验：
  - evidence 的 `source` 必须存在于 lineage；
  - evidence scope 不得超出 source analysis 的 scope；
  - **AI claim 不能被再次当作独立 deterministic evidence 引用**。

**验收**：

- 测试：跨路线共享输入的 evidence 被标记为同一 `dependency_group`；
- 测试：scope 越界被拒绝；
- 测试：AI claim 互引被拒绝。

---

### T5. Strict privacy 收口 + 用户 consent UX

**目标**：把 privacy 从「默认在 mock 之外启用」变成**明确的、可审计的**行为。

**做法**：

- `sclet_ai_call()` 的 privacy 行为明确化：
  - 真实 provider 调用：默认 `standard` 且**真实执行**裁剪；
  - strict 模式遇到未分类字段**拒绝**（现有 `sclet_ai_payload_scan` 已实现，需保证在真实调用路径上生效）；
  - `privacy_consent = TRUE` 的语义写清楚（哪些字段会因 consent 放行、哪些是**永久硬拒绝**：token/key/path/matrix）。
- 增加**审计日志字段**：每次出站 payload 的 `payload_fingerprint`、`redacted_fields` 必须可从结果 `metadata$payload_policy` 读到（已部分实现，需确认真实路径也写入）。
- consent UX：当裁剪导致 `requires_user_consent = TRUE` 时，返回可读提示，**不要静默丢字段**。

**验收**：

- 测试：真实调用路径（mock provider，但 `enforce_privacy = TRUE`）不含 secret / matrix / raw ID；
- 测试：strict 未分类字段被拒绝；
- 测试：`payload_fingerprint` 稳定且可审计；
- 测试：consent 不能放行 token / path / matrix（永久硬拒绝）。

---

### T6. 文档与 Spec 同步

**做法**：

- 手写 `man/*.Rd`（不要 `make rd`）；
- 手工在 `NAMESPACE` 加 `export(...)`（注意**不要重复**，之前出过重复行）；
- 手工在 `DESCRIPTION` 的 `Collate:` 插入新 R 文件（顺序：被依赖的先加载）；
- `NEWS.md` 顶部加一条本阶段条目（中文可，但注意 NEWS.md 不是 R 源码，中文没问题）；
- 更新 `.dev/ai-advanced-analysis-spec.md`：
  - 状态行、`重要状态说明`、`还需要建设`、第一批/下一批任务列表；
  - **只写真实已实现的能力**，proposed API 必须标注未实现。

---

## 4. 验收清单（交付前必须全部满足）

```bash
# 1. 纯 ASCII
LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-*.R && echo "FAIL: non-ASCII"

# 2. 无空白错误
git diff --check

# 3. 聚焦 AI 测试全绿（在线测试应 skip）
Rscript -e 'devtools::test(filter = "ai-", reporter = "progress")'

# 4. 完整包检查干净
make check
```

目标：

```text
[ FAIL 0 | WARN 0 | SKIP 1 | PASS >= 223 ]     # SKIP 1 = 在线测试
0 errors | 0 warnings | 0 notes                  # make check
```

**交付物**：

- [ ] 新增/修改的 `R/*.R`
- [ ] 新增/修改的 `tests/testthat/test-ai-*.R`
- [ ] 手写的 `man/*.Rd`
- [ ] `NAMESPACE` / `DESCRIPTION` / `NEWS.md` 更新（无重复、无无关改动）
- [ ] `.dev/ai-advanced-analysis-spec.md` 状态同步
- [ ] 最终报告：**贴出上面 4 条命令的真实输出**，并逐条说明 T1–T6 的完成情况与**未完成项及原因**

---

## 5. 已知陷阱（踩过，别再踩）

1. **`make rd` 会破坏 NAMESPACE** —— 手工维护。
2. **R 源码非 ASCII 会 WARNING** —— 用 `\uXXXX`。
3. **`R/ai-*` 文件里用 `capture.output` 必须写 `utils::capture.output`**，否则 check 出 NOTE。
4. **`sclet_set_analysis_state()` 只接受 `sclet_state_types()` 里的 type** —— 新增 type 必须同时改 `R/state.R` 的 `sclet_state_types()`（`ai_evidence` 已加）。
5. **ledger 里同一分析会出现两份**（accessor record + state record），做比较/查 source 时必须按 `id` 去重；`state_records` 是嵌套 list，需要 `unlist(..., recursive = FALSE, use.names = FALSE)` 展平。
6. **harmony 未安装** —— 不要写依赖真装的测试；不要把它变成硬依赖。
7. **`%||%` 是本包内部工具**，直接用即可。
8. **非交互式会话 `confirm = "ask"` 不执行**，只回 dry-run —— 这是刻意设计，不要"修好"它。
9. **不要把 `recommendation` 填成"最佳路线"** —— 指标冲突时必须留给用户决定。
10. **basilisk 类函数用 mock 测试**，不要真的建 Python 环境。

---

## 6. 完成后

把最终报告（命令输出 + T1–T6 逐条状态 + 未完成项）留在仓库或回复里。原 agent 会按以下角度复核：

- 是否**真的**执行了路线，还是只登记了摘要；
- 是否有**未经确认就猜 batch/condition** 的路径；
- 是否存在**单分数自动选优**；
- evidence / scope / fingerprint 校验是否可被绕过；
- privacy 是否在**真实调用路径**生效（而不只是 mock 之外）；
- 测试是否**真的能失败**（不是恒真断言）；
- `make check` 是否 `0/0/0`。
