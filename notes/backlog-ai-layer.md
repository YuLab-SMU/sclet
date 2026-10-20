# AI 层小问题台账（2026-10-19 起）

> 用途：**记录，不追着修**。开工先读 `.dev/ai-product-roadmap.md` 定里程碑，
> 再看这里有没有与该里程碑相关、且已经阻塞结论的条目。其余一律往后放。
>
> 分级只按一个问题回答：**它会不会污染我们即将说出口的结论？**
>
> - `claim` —— 会。当场修，并配 pin 测试（"a fix without a pinned test is a claim"）。
> - `report` —— 不污染数据，但会让报告/账本说出不实或不完整的话。攒到报告层开工时一起修。
> - `noise` —— 只影响阅读体验或开发效率。顺手修，不专门开一轮。
>
> 每条都要有「回访触发条件」。没有触发条件的条目不是台账，是愿望。

## 已修（2026-10-19 真实数据跑出来的六个，均为 `claim` 级）

PBMC 10k（11,769 × 33,538）外部 agent 跑分轮计划时发现，细节见 `NEWS.md` 当日条目
和 `tests/testthat/test-ai-realdata-provenance-regressions.R`：

1. marker 证据节点记录写死的阈值（0.05 / 0.1）而非本次运行的阈值 → 现已读参数。
2. `set_dimred()` / `RunUMAP()` 把 `.dimred` 永久留在对象上（11,769×30 副本）→ 现已不留。
3. `run_enrichment` 缺 `db` 时静默用 GO → 现在 `db_missing` / `db_invalid` 拒绝。
4. 账本不暴露 HVG 选择、执行记录不含每步参数与耗时 → 现在有 `features` 与 `params$steps`。
5. KEGG 空结果（`n_enriched = 0`、`result` 为 NULL）仍报 `completed` → 现在 `enrichment_no_result` 失败。
6. louvain 聚类无种子、不可复现（同一计划 13 群 vs 12 群）→ 现在 `seed` 参数并记录。

## 已修（2026-10-20 M3：报告层）

7. **B1 冲突误报**：`analysis_story$conflicts` 从第二轮起必然报 `ai_execution`。现在追加式日志类型
   （`ai_execution` / `ai_evidence` / `ai_design_confirmation`）不计入冲突；计算型状态（detest、
   clustering、rare_cells…）仍然报告。pin：`test-ai-report-layer.R` 第 1 条。
8. **B2/B3 marker 证据的假字段**：删掉 `top_up_ranks`（构造上恒为 1:10）、`top_up_fractions`
   （恒为 10 个 NA）、两个无消费者的不透明整数编码，以及 `n_significant`（表已被调用方阈值
   预过滤，再数一遍等于复述过滤；实测 12,559/12,572）。`n_top_up` 改成真实计数。
   pin：`test-ai-report-layer.R` 第 2 条。
9. **B4 detest 记录缺参数**：`RunDEtest()` 现在把 `min.pct` / `logfc.threshold`（以及条件路径的
   `min_cells`）写进 detest 记录的 `params`，单独读 `get_detest()` 就能复跑。pin：同文件第 3 条。
10. **B14 报告渲染层**：新增导出 `RenderAIReport()`（`R/ai-report.R`）。从账本渲染 Markdown，
    并**自检**：引用必须可解析、每条记录等级必须等于引擎重算值，不一致就打印在「复验」区而不是
    抹平。不调用模型、不含原始值。M2 报告现在可以由它从账本重现（`m2-verification-kang.md` V12/V13）。
11. **B10 只读步骤的痕迹**：`check_qc` 这类只读动作仍不写分析记录（只读契约要求不 mutate），
    但执行记录的 `params$steps` 有它的一行，`RenderAIReport()` 的「执行溯源」区会列出来。
    所以「计划里跑过 QC」现在可查。pin：`test-ai-report-layer.R` 第 5 条。
12. **结果自带的 warnings 曾丢失**：`RecordAIResult()` 的 artifacts 里没有 `warnings`，报告渲染
    时看不到运行提出的告诫（M2 的 3 条 warning 只活在内存里）。现在随结果入账。
    pin：`test-ai-report-layer.R`（warnings 出现在渲染结果里）+ M2 复验 V13。

13. **派生 marker 表上的富集丢了 parent（claim 级）**：KEGG 富集跑在 bitr 转换后的派生 detest 记录上，
    那条记录自己没有证据节点，`sclet_ai_enrichment_parent_nodes()` 只按直接 source 查 → KEGG 节点
    没有 parents → `check_claim_ceiling()` 把「markers + GO + KEGG」算成**两条独立证据线**，
    而它们本是同一份 marker 表的两种注释。现在会跟一层 `inputs$parent_detest`。
    真实数据复验：M1 全量证据的 `n_independent_groups` 从 2 变回 1。
    pin：`test-ai-report-layer.R` 第 6 组；M1 复验 V10。

## 未修

### B1 `analysis_story$conflicts` 从第二轮起必然报警 · `report`
`"multiple ai_execution records exist; only the active one (if any) is in active_states"`
对任何多轮分析都成立——`ai_execution` 本来就是追加式日志。总在响的冲突标志会让读者
学会忽略冲突字段。建议：把 `ai_execution`（以及 `ai_evidence`）排除在"冲突"之外，
或改写成 `history` 语义。
- 证据：`/tmp/sclet-agentrun/31-deep.R` 在 r2 上打印 1 条 conflict。
- 回访触发：报告层开始渲染 conflicts 时。

### B2 marker 证据里几个字段名不副实 · `report`
`top_up_ranks` 恒为 `1:10`（不是名次），`top_up_fractions` 恒为 10 个 `NA`，
`annotation_method_code` / `annotation_scope_code` 是不透明整数。都是"承诺了但没测量"的形状。
- 回访触发：证据要进入报告正文时。

### B3 `n_significant` 几乎等于 `n_rows` · `report`
`FindAllMarkers` 的返回表已经是按 `p_val_adj` 过滤过的，所以再按 0.05 数一遍
几乎必然全中（实测 12,559 / 12,572）。这个数字读起来像"99.9% 显著"，实际是
"表已经过滤过"。建议改名（`n_rows_reported`）或直接删掉。
- 回访触发：同上。

### B4 detest 记录的 `params` 仍为空 · `report`
阈值现在能在两处看到（证据节点 `values`、执行记录 `params$steps`），但 detest 记录本身
没有，单独读 `get_detest()` 的人看不到阈值。`RunDEtest()` 加可选参数即可，但会动公开 API。
- 回访触发：有人只凭 `get_detest()` 复跑时。

### B5 没有 symbol→ENTREZID 的 action，KEGG 在计划层面不可达 · `product`
今天跑 KEGG 必须在校验层之外先 `bitr`，再把转换后的表注册成一个派生的 detest 记录，
计划才跑得通。这既是缺口也是真实用户会踩的坑。可选：加一个 `convert_gene_ids` action，
或在 `run_enrichment` 的 KEGG 分支里做转换并记录。
- 回访触发：下一次出现 KEGG/Reactome 之类的目标。

### B6 账本没有构建身份 · `dev`
版本号长期停在 1.0.2，同一版本号下 `~/Rlib` 可能是几天前的构建。今天实际踩到：
先跑了两天前的 build，靠 schema 差异才发现。建议 `dataset` 或顶层带
`sclet_version` + 构建时的 git sha。
- 回访触发：下一次"这不是已经修了吗"的时刻。

### B7 `sclet-analysis` 的表格断行 · `noise`
`inst/skills/sclet-analysis/SKILL.md` 99–102 行：四行表格行落在两段正文之后，
渲染不成表格。
- 回访触发：下次改这个文件时顺手修。

### B8 旧对象上的 `.dimred` 没有迁移 · `noise`
新代码不再产生 `.dimred`，但旧 `saveRDS` 的对象里还留着。载入旧 fixture 时可能看到
`reducedDimNames` 多一项。
- 回访触发：复用到 2026-10-19 之前保存的对象时。

### B9 载入期噪声 · `noise`
`library(sclet)` 之后首次用到 enrichment 会打印 `Loading required namespace: enrichit`，
以及 `Found more than one class "package_version"`（SeuratObject vs alabaster.base）。
不影响结果，但污染脚本日志，也让人以为是错误。
- 回访触发：做批量运行时。

### B10 `check_qc` 不留分析记录 · `report`
只读动作返回一个值，账本里没有对应记录；今天补的执行记录 `params$steps` 是它唯一的痕迹。
"计划里跑过 QC"这件事在账本里无法独立查到。
- 回访触发：报告需要写"做过哪些检查"时。

### B11 `observed` 表达不了"对象结构事实" · `report`
写报告时发现：一条不带 `evidence_refs` 的结构事实（"对象有 11,769 个细胞"）会被审计
降级为 `hypothesis`——`observed` 要求"事实在对象上"，但审计只认可解析的证据 id，而
"对象形状"没有对应证据节点。M1 的做法是把这类事实放进报告的「方法与参数」正文，不当作
finding。若要让结构事实以 `observed` 入账，需要一条对象状态证据节点，或允许 `observed`
绑定状态记录 id。
- 证据：`/tmp/sclet-agentrun/100-claim-probe.R`，`observed` + 无引用 → `derived=hypothesis`，
  `rule=no_prose_reference`。
- M3 后：报告层不再需要它——结构事实进「执行溯源 / 方法与参数」区（描述性，不评级），
  只有带证据的结论才进 findings。真正的 `observed` 入账仍需设计（对象状态证据节点）。
- 回访触发：有人确实需要把"对象有 N 个细胞"作为 observed finding 入账时。

### B12 审计与写入对"未知证据 id"行为不一致 · `report`
`sclet_ai_audit_result_claims()` 把未知 id 静默降级为 `hypothesis`（理由写
`unknown or incomplete evidence id(s)`），而 `RecordAIResult()` 在同一情形直接拒绝。
只调审计的调用方会以为"等级被校正"，实际是引用写错了。
- 证据：`100-claim-probe.R`（写错 `..._kegg_entrez` → 降为 hypothesis）；
  `102-record-probe.R`（同形状经 `RecordAIResult` → `REFUSED: evidence ref not valid`）。
- 回访触发：下一次有人只用审计 API 做校验时。

### B13 账本默认丢 artifacts，读回自己的结果会静默为空 · `noise`（M3 后：靠约定避开）
`RenderAIReport()` 内部固定 `include_artifacts = TRUE`；任何读回自己写的 AI 结果的代码都必须显式传。
`GetAnalysisLedger(detail = "full")` 默认 `include_artifacts = FALSE`，所以刚记入的
AI 报告的 `artifacts$findings` 读回来是 `NULL`，`length(...)` 为 0 且不报错（M1 生成脚本
第一次就打印了 `recorded findings: 0`）。建议：要么文档写明，要么在 detail = "full" 下
默认带上 artifacts 的元数据。
- 回访触发：任何"读回来怎么是空的"时刻。

### B14 没有报告渲染层 · `product`
report contract（`.dev/ai-product-roadmap.md` §3.3）里的 Analysis Report 目前只以
`sclet_ai_result` 记录存在，包内没有把它渲染成可读产物的函数；M1 的 markdown 由驱动脚本
手工拼装。若报告要作为产品出口，需要 `RenderAIReport()` 之类的一层（输入 = 记录 id +
对象，输出 = markdown/html），否则每个驱动者都会各写一份渲染器。
- 回访触发：M2 收尾前决定。

### B15 两个证据访问器形状不一致 · `noise`
`sclet_ai_evidence_get_all()` 返回裸节点（`$values` 在顶层），账本的
`state_records$ai_evidence[[id]]` 包在 `$summary` 里。取错会静默得到 `NULL`——M2 生成脚本
因此写出过一条空数字的 finding，被复验拦下。建议：给裸节点补一个 `$summary` 别名，或在
`sclet_ai_evidence_get_all()` 文档里写明形状。
- 证据：`/tmp/sclet-agentrun/220-m2-report.R` 的修复注释；NEWS 当日条目末段。
- M3 后：`RenderAIReport()` 与复验脚本统一用账本形式（`state_records$ai_evidence[[id]]$summary$values`），
  并在 `sclet_ai_evidence_get_all()` 的注释里写明它返回的是**summary 扁平化**后的节点。
- 回访触发：下次有人取证据节点的 values 时（大概率还会踩）。

### B16 scVI 全量对象补跑（43,095 细胞） · `product`
M2 的 scVI 只在 9,000 细胞子集上跑（2,878s 含首次 env 构建）。env 已建好，全量对象随时可跑
（CPU 上数小时量级）。跑完把 scVI 行并入 `notes/m2-report-kang.md` 的路线表，并重跑
`compare_integration_routes`（会得到 4 条路线的新证据节点）。
- 证据：`/tmp/sclet-agentrun/211-m2-scvi.R`（子集脚本，改 keep 的每组上限即可全量）。
- 回访触发：路线比较表需要第三条路线时，或有人问"scVI 全量是多少"。

## 维护约定

- 新增条目只写三件事：现象 + 证据位置 + 回访触发条件。
- 修掉的条目移到「已修」并写明 pin 测试文件与断言，不删原文。
- 一次里程碑收尾时过一遍台账，只做与下个里程碑相关的条目。
