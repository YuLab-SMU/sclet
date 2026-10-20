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

## 维护约定

- 新增条目只写三件事：现象 + 证据位置 + 回访触发条件。
- 修掉的条目移到「已修」并写明 pin 测试文件与断言，不删原文。
- 一次里程碑收尾时过一遍台账，只做与下个里程碑相关的条目。
