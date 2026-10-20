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

14. **M4 单一报告来源**：`notes/m1-report-pbmc10k.md`、`notes/m2-report-kang.md` 改为
    `RenderAIReport()` 输出；渲染器新增「Recorded analyses」（含路线比较表）与
    「Limitations stated in the findings」两区。手写叙述中无法从账本重建的部分：
    scVI 子集尝试（另一个对象）→ 归 B16，报告不再提；其余全部由记录渲染。

15. **M5 annotation 闭环**：PBMC 10k 的 13 个群完成参考注释（MonacoImmuneData / label.main，
    9.5 s）；参考不是引擎猜的，而是走**澄清通道**记录成 `user_decision` 后由下一轮计划读取
    （`ResolveAIClarifications` 全链路实测）。markers + annotation 触发 2 条独立证据线，
    上限仍 `associated`（最弱节点）。顺带修掉 annotation 证据节点里同类假字段：
    `top_label_ranks`（构造上恒为 1:k）与两个不透明整数编码。
    报告：`notes/pbmc10k-annotation-report.md`（`RenderAIReport()` 输出）+ 复验 9/9。

16. **M6 milo 闭环（差异丰度）**：Kang 上 25 个样本单元（donor × condition）、3,414 个邻域，
    stim vs ctrl 有 1,821 个 FDR<0.05；**置换标签对照 0/3,414（中位 p 0.76）**，所以信号不是
    标签顺序造成的。设计 `~ 0 + condition + batch`（未配对）作为限制写进报告。
    顺带修掉两个“从没被执行过所以没人发现”的缺陷：
    - `run_milo` 的错写形式（技能与错误信息都教 `"stim - ctrl"`）被原样送进
      `miloR::testNhoods()`，而它按设计矩阵系数名解析 → `object 'stim' not found`。
      现在 `RunMilo()` 在知道设计的地方把 level 对比翻译成系数名（`conditionstim - conditionctrl`），
      记录里同时保留调用方写法与模型写法；pin：`test-ai-report-layer.R`。
    - action handler 假设 `RunMilo()` 返回 SCE，实际返回 `list(sce=, milo=, da_results=)` →
      每次执行都在最贵的那一步之后死于 `no metadata method for class list`。
      覆盖率表里 `run_milo` 一直是 `none`，这正是它没被发现的原因。

17. **M7 doublet / rare 闭环（并回答了 M5 留下的问题）**：PBMC 10k 上 scDblFinder 判出
    863/11,769 双细胞（7.33%，467 s）。两个具体结论：
    - **cluster 12 是双细胞群（291/296，98.3%）**，它在 M5 里被 Monaco 标成 Monocytes——
      这是对先前组成行的**修正**，M7 报告里明说，M5 的记录不改（记录是历史）。
    - **cluster 9 的 doublet 是 0/719**，所以 M5 里"过渡态或双细胞"的假设排除了双细胞那一半。
    - 稀有群：2,350 个密度极大值解散后 **0 个群体**，按设计不注册证据节点，报告把它写成
      hypothesis 而不是结论。
    报告 `notes/pbmc10k-doublet-report.md`，复验 8/8（含"零稀有群 ⇒ 无证据节点"这条规则）。

18. **M8 cci 闭环（引擎由用户选定）**：CellChat 在 12 个发送群 / 13 个接收群之间推断出 733 对
    显著互作（113 s），引擎经澄清通道记录为 `user_decision`（第二次使用该通道）。
    报告写明两件事：① 本机只装了 CellChat，所以这个“用户选择”也受环境限制；
    ② **cluster 12（M7 里 98.3% 双细胞）参与了互作推断**，涉及它的互作含污染。
    证据节点不含配体/受体/通路名称（`method_code` 是有意的匿名编码，不是死字段——已核对注释）。
    报告 `notes/pbmc10k-cci-report.md`，复验 9/9。
    复验脚本自己也踩过一个坑：用 `IL` 做泄漏检查会命中 `probabILity`，改成词边界后才有意义。

19. **M9 轨迹链 + inspect_ledger**：`inspect_ledger` 已执行（最后一个 `none` 动作清掉）。
    真实数据上 slingshot 得到 1 条 lineage（起点 cluster_1 = Pachytene，5 群 300 细胞）。
    **velocity 在真实数据上被数据挡住**：velociraptor 需要 `counts` assay，而本仓库所有
    velocity 导出（hermann sperm、regvelo smoke）只有 spliced/unspliced/logcounts →
    计划层现在直接拒绝（`counts_assay_missing`），不再等到 Python 里崩。
    速度 + fate 在**合成 fixture** 上跑通（scvelo deterministic + CellRank 43 s），但 fate
    **退化为单一终态**（1 个 terminal state、吸收概率全 1）——链路通不等于结果有意义，
    报告如实标注。附带成本：cellrank 的 basilisk env 按包版本号存路径，1.0.1 → 1.0.2
    升版后要重建（一次 ~10 分钟）。

20. **收尾这一轮修的四个 report 级问题**：
    - CCI 证据的 `n_significant_pairs` 在 CellChat 上恒等于行数（表已被引擎按 p<0.05 过滤，
      实测 733/733、max p = 0.04）。新增 `p_value_max`，让读者能分辨“表已过滤”与“真的都显著”，
      字段语义因此成立。
    - 审计对未知证据 id 以前只把等级降成 `hypothesis`，而 `RecordAIResult()` 会直接拒绝——
      只审计不写入的调用方看到的是“安静降级”。现在 entry 带 `unresolved_evidence_refs` 标志，
      原因可见；pin 在 `test-ai-report-layer.R`。
    - B13/B15 写进 `GetAnalysisLedger()` 文档与 `sclet-driver` 技能。
    - B11 写进 `sclet-claims` 技能：结构事实属于出处/方法，不属于 findings（无引用的 finding
      必然降为 hypothesis）。

21. **M10 CNV / malignant 闭环（并复现了用户反馈的两个结构效应）**：GSE144735 KUL3（1,400 细胞，
    300 参考），inferCNV HMM on 跑 **61.8 min**；两份四态调用（纯双线 / 加免疫竞争签名）。
    CRC 上的数字与 GSE150290 的报告**同构**：
    - 签名阈值取了全池中位数（-0.00893），签名阳性 **700/1400 = 恰好一半**——这是定义，不是生物学；
    - CNV 阈值 0.183（参考 P99+3×MAD），CNV+ 只有 **7/1400**；
    - 加 identity 竞争后 `malignant` 从 **4 → 0**（降级前 supported 4 保留在
      `<name>_call_pre_demote`），`identity_conflict = 1`。
    顺带修了三件事：CNV 证据节点记录 `parameter_hmm`/`parameter_cutoff`/`parameter_n_genes`；
    malignant 节点记录 `signature_cut_*`/`cnv_reference_threshold`/`n_cnv_positive`/
    `n_pre_demote_supported`；`RecordAIEvidence(replace = TRUE)` 允许在**不重跑分析**的前提下
    重记派生节点（inferCNV 一小时一次，旧节点无法更新是缺陷）。
    报告 `notes/crc-malignant-report.md`，复验 10/10；对 Issue #26 的核查见
    `notes/malignant-gse150290-analysis-2026-10-20.md`。

22. **M11 真实数据的 速度-轨迹-命运 链（补上最后一个数据缺口）**：用 scvelo 的 pancreas 数据
    （Bastidas-Ponce 2019，3,696 细胞，带 counts/spliced/unspliced）——真实数据，不是 fixture。
    ``run_velocity`` deterministic 45 s；slingshot 从 Ngn3 low EP 出发得 **5 条 lineage**；
    CellRank 9 min 得 **3 个终态**。
    顺带发现并修掉一个 report 级缺陷：证据里的 ``terminal_state_fractions`` 是**终态细胞之间**
    的占比（90/3,696，按 1/k 构造，实测 0.333/0.333/0.333），读起来却像"每个命运去多少细胞"。
    现在节点同时记录 ``n_terminal_state_cells``、``destination_fractions``（全细胞吸收概率
    argmax，实测 **0.107/0.599/0.294**）与 ``fate_dimension_mean_absorption``，
    pin 在 ``test-ai-fate.R``。
    报告 ``notes/pancreas-velocity-fate-report.md``，复验 10/10。

23. **重跑后证据节点静默保留旧数字（2026-10-21 自查发现，未在台账里）**：
    节点 id 从分析名派生（`ev:doublet_doublet_finder`），所以同一动作跑两次会撞同一个 id；
    `RecordAIEvidence()` 拒绝重复，而**七个 recorder 全都**用
    `tryCatch(error = function(e) object)` 把它吞掉。实测：colData 里 0 个 doublet，
    节点仍写 150；scDblFinder 是随机的，真的第二次跑会得到不同数字（fixture 上 4 → 2）。
    这是"引用了过时数字"的缺陷——比没有节点更糟，因为它可引用、会渲染、看起来像本次结果。
    修法是一处共享函数 `sclet_ai_record_evidence_for_run()`，把三种诚实的结局分清：
    - 数字相同 → no-op（幂等，不报错）；
    - 数字不同且未要求 replace → **报 `evidence_id_conflict`**，说明两次运行不一致；
    - `replace = TRUE` → 刷新节点。
    7 个 recorder（doublet / cnv / milo / velocity / trajectory / fate / enrichment /
    communication / routes / malignant）全部改走这条路径，并在 `test-ai-evidence-rerun.R`
    里逐函数断言"必须路由到它"，所以以后新增 recorder 漏掉会被测出来。
    顺带修掉一个**把 bug 当预期固定下来的测试**：`test-ai-interpretation-gaps.R` 原来断言
    "第二次真实运行必须降级为 no-op"，现已改成断言拒绝 + replace 刷新。

## 未修

### B1 `analysis_story$conflicts` 从第二轮起必然报警 · `report`（2026-10-20 关闭：M3）
`"multiple ai_execution records exist; only the active one (if any) is in active_states"`
对任何多轮分析都成立——`ai_execution` 本来就是追加式日志。总在响的冲突标志会让读者
学会忽略冲突字段。建议：把 `ai_execution`（以及 `ai_evidence`）排除在"冲突"之外，
或改写成 `history` 语义。
- 证据：`/tmp/sclet-agentrun/31-deep.R` 在 r2 上打印 1 条 conflict。
- 回访触发：报告层开始渲染 conflicts 时。

### B2 marker 证据里几个字段名不副实 · `report`（2026-10-20 关闭：M3）
`top_up_ranks` 恒为 `1:10`（不是名次），`top_up_fractions` 恒为 10 个 `NA`，
`annotation_method_code` / `annotation_scope_code` 是不透明整数。都是"承诺了但没测量"的形状。
- 回访触发：证据要进入报告正文时。

### B3 `n_significant` 几乎等于 `n_rows` · `report`（2026-10-20 关闭：M3）
`FindAllMarkers` 的返回表已经是按 `p_val_adj` 过滤过的，所以再按 0.05 数一遍
几乎必然全中（实测 12,559 / 12,572）。这个数字读起来像"99.9% 显著"，实际是
"表已经过滤过"。建议改名（`n_rows_reported`）或直接删掉。
- 回访触发：同上。

### B4 detest 记录的 `params` 仍为空 · `report`（2026-10-20 关闭：M3）
阈值现在能在两处看到（证据节点 `values`、执行记录 `params$steps`），但 detest 记录本身
没有，单独读 `get_detest()` 的人看不到阈值。`RunDEtest()` 加可选参数即可，但会动公开 API。
- 回访触发：有人只凭 `get_detest()` 复跑时。

### B5 没有 symbol→ENTREZID 的 action，KEGG 在计划层面不可达 · `product`（2026-10-20 关闭：KEGG+SYMBOL 现在验证期就拒绝，并给出两条出路）
今天跑 KEGG 必须在校验层之外先 `bitr`，再把转换后的表注册成一个派生的 detest 记录，
计划才跑得通。这既是缺口也是真实用户会踩的坑。可选：加一个 `convert_gene_ids` action，
或在 `run_enrichment` 的 KEGG 分支里做转换并记录。
- 回访触发：下一次出现 KEGG/Reactome 之类的目标。

### B6 账本没有构建身份 · `dev`（2026-10-20 关闭：账本带 `sclet_version` 与安装时间）
版本号长期停在 1.0.2，同一版本号下 `~/Rlib` 可能是几天前的构建。今天实际踩到：
先跑了两天前的 build，靠 schema 差异才发现。建议 `dataset` 或顶层带
`sclet_version` + 构建时的 git sha。
- 回访触发：下一次"这不是已经修了吗"的时刻。

### B7 `sclet-analysis` 的表格断行 · `noise`（2026-10-20 关闭：补回表头与空行）
`inst/skills/sclet-analysis/SKILL.md` 99–102 行：四行表格行落在两段正文之后，
渲染不成表格。
- 回访触发：下次改这个文件时顺手修。

### B8 旧对象上的 `.dimred` 没有迁移 · `noise`（won't fix：旧对象重跑即可，不写迁移代码）
新代码不再产生 `.dimred`，但旧 `saveRDS` 的对象里还留着。载入旧 fixture 时可能看到
`reducedDimNames` 多一项。
- 回访触发：复用到 2026-10-19 之前保存的对象时。

### B9 载入期噪声 · `noise`（won't fix：来自 SeuratObject/alabaster 的 S4 类重名，不在本包控制内）
`library(sclet)` 之后首次用到 enrichment 会打印 `Loading required namespace: enrichit`，
以及 `Found more than one class "package_version"`（SeuratObject vs alabaster.base）。
不影响结果，但污染脚本日志，也让人以为是错误。
- 回访触发：做批量运行时。

### B10 `check_qc` 不留分析记录 · `report`（2026-10-20 关闭：M3 渲染器的「Executed steps」）
只读动作返回一个值，账本里没有对应记录；今天补的执行记录 `params$steps` 是它唯一的痕迹。
"计划里跑过 QC"这件事在账本里无法独立查到。
- 回访触发：报告需要写"做过哪些检查"时。

### B11 `observed` 表达不了"对象结构事实" · `report`（2026-10-20 关闭：写进 skill，不改阶梯）
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

### B12 审计与写入对"未知证据 id"行为不一致 · `report`（2026-10-20 关闭：审计新增 `unresolved_evidence_refs` 标志）
`sclet_ai_audit_result_claims()` 把未知 id 静默降级为 `hypothesis`（理由写
`unknown or incomplete evidence id(s)`），而 `RecordAIResult()` 在同一情形直接拒绝。
只调审计的调用方会以为"等级被校正"，实际是引用写错了。
- 证据：`100-claim-probe.R`（写错 `..._kegg_entrez` → 降为 hypothesis）；
  `102-record-probe.R`（同形状经 `RecordAIResult` → `REFUSED: evidence ref not valid`）。
- 回访触发：下一次有人只用审计 API 做校验时。

### B13 账本默认丢 artifacts，读回自己的结果会静默为空 · `noise`（2026-10-20 关闭：写进 GetAnalysisLedger 文档 + driver skill）
`RenderAIReport()` 内部固定 `include_artifacts = TRUE`；任何读回自己写的 AI 结果的代码都必须显式传。
`GetAnalysisLedger(detail = "full")` 默认 `include_artifacts = FALSE`，所以刚记入的
AI 报告的 `artifacts$findings` 读回来是 `NULL`，`length(...)` 为 0 且不报错（M1 生成脚本
第一次就打印了 `recorded findings: 0`）。建议：要么文档写明，要么在 detail = "full" 下
默认带上 artifacts 的元数据。
- 回访触发：任何"读回来怎么是空的"时刻。

### B14 没有报告渲染层 · `product`（2026-10-20 关闭：`RenderAIReport()`）
report contract（`.dev/ai-product-roadmap.md` §3.3）里的 Analysis Report 目前只以
`sclet_ai_result` 记录存在，包内没有把它渲染成可读产物的函数；M1 的 markdown 由驱动脚本
手工拼装。若报告要作为产品出口，需要 `RenderAIReport()` 之类的一层（输入 = 记录 id +
对象，输出 = markdown/html），否则每个驱动者都会各写一份渲染器。
- 回访触发：M2 收尾前决定。

### B15 两个证据访问器形状不一致 · `noise`（2026-10-20 关闭：写进 driver skill）
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
- 命令：在 `m2-final.rds` 上 `RunIntegration(method = "scVI", batch = "batch", name = "scvi")`
  再 `sclet_ai_record_route_metrics()`；env 已按 1.0.2 建好（M9 顺带重建了 cellrank env）。
- 回访触发：需要一个覆盖全部 43,095 细胞的第三条路线时。

## 维护约定

- 新增条目只写三件事：现象 + 证据位置 + 回访触发条件。
- 修掉的条目移到「已修」并写明 pin 测试文件与断言，不删原文。
- 一次里程碑收尾时过一遍台账，只做与下个里程碑相关的条目。
