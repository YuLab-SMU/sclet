# 对 Issue #26 最新回复的核查（2026-10-20）

> 对象：`YuLab-SMU/ProjectYulab#26` 中 `YoupengFan` 于 2026-10-20 03:59 的回复
> （GSE150290 上的对照实验：迭代基线 + 双线四态 + DNA 仲裁 + cutoff 敏感度）。
> 结论：**这是一次高质量对照，指出的三个问题在 sclet 里都成立**，其中两个是默认参数的
> 语义问题、一个是证据层缺字段导致无法复算。本文逐条对码核查，并记下要改什么。

## 一、三条主张，逐条对码

### 1. 「签名+ 恰为全池一半」——成立，而且是我们默认参数的直接后果

`R/phenotype-evidence.R` 的 `transcriptomic_cut` 默认取**全池分数中位数**：

```r
cut_used <- if (is.null(transcriptomic_cut)) {
    list(value = stats::median(sig_score, na.rm = TRUE), source = "median")
} else {
    list(value = transcriptomic_cut, source = "caller")
}
sig_positive <- !is.na(sig_score) & sig_score >= cut_used$value
```

所以 44,011 个细胞里 22,006 个「签名+」不是数据特征，是**定义**（`>= 中位数`）。
用户把它写成「语义是"种子签名比一半细胞激活得强"的相对陈述」，这个读法完全正确，
而且我们的代码里 `source = "median"` 就是这个意思——但它只活在分析层，**没有进证据节点**，
所以从账本里看不出这 22,006 是定义而非测量。

同理，迭代法 24,390（55.4%）也来自 k-means 的 k=2 相对二分。我们之前已经为这个形状
加了 `split_degenerate` 警告，但**恶性率接近 50% 这个事实本身**没有被当成异常信号。

### 2. 「identity 竞争检验把全部肯定结论降级」——成立，而且是数据集级的

`R/phenotype-evidence.R`：

```r
demoting <- names(Filter(function(x) isTRUE(x$demote), competing_report))
if (length(demoting) > 0) {
    call[call %in% c("supported", "probable")] <- "undetermined"
}
```

一个数据集级的分离度（用户实测 +0.199，阈值 0.1）就把**所有** supported/probable 改写为
undetermined。用户看到的「29 个双阳 → 全部 undetermined → malignant 0」正是这条语句。
后果是个结构性的：只要竞争签名在这个数据集上分离超过阈值，`malignant` 这个状态就**不可达**，
四态塌成两态。这不是生物学结论，是设计后果。

### 3. 「CNV 线在 DNA 真值上零命中」——成立，但要分清是阈值还是队列

用户的金标准：Kim et al. 2022 Supplement Table 6 的 201 个热点突变细胞，**没有一个 CNV+**
（全池 29 个 CNV+ 都不在其中）。我们的 CNV 阈值是

```r
base <- stats::quantile(ref_burden, reference_quantile, names = FALSE)   # 默认 P99
candidate <- base + 3 * stats::mad(ref_burden)
```

即「参考细胞 burden 的 P99 + 3×MAD」。在 10,853 个控制细胞上，这是一个很严的阈值——
44,011 个细胞只有 29 个越过。所以两种解释同时成立：**阈值确实偏保守**，而且**这个队列的
恶性可能以点突变为主、CNV 稀少**（用户自己的判断）。但它们指向不同的修法：
前者要报告灵敏度的不确定性，后者要换证据线。

## 二、一个必须先说清的前提：他们跑的不是 HMM 版

用户明确写了 `hmm = FALSE`（32 GB 内存 OOM），全量改用无 HMM 的窗口平滑版。
我们文档里描述的 CNV 路径默认 `hmm = TRUE`，burden 来自 HMM（i6）profile。
**两者不是同一个对象。**

而我们的证据节点当时**根本不记录 `hmm`**——所以他们的运行与 HMM 运行在账本里无法区分。
这是我们的缺陷，不是他们的表述问题。（已在本次修掉，见下。）

其余两条实现说明都站得住：`cutoff = 0.5` 预子集 980 基因（对照 `cutoff = 1` 的 524 基因，
结论只差 1 个细胞）、以及 44 个边界细胞随阈值翻动——这正是「阈值附近的细胞不该被当成稳定结论」
的证据。

## 三、这次改了什么（对应上面的三条）

1. **证据节点记录判决的「工作点」**（`R/ai-execution.R`）：
   - CNV 节点增加 `parameter_hmm`、`parameter_cutoff`、`parameter_n_genes`；
   - malignant 节点增加 `signature_cut_value`、`signature_cut_is_median`、
     `n_signature_positive`、`cnv_reference_threshold`、`n_cnv_positive`、
     `n_pre_demote_supported`。
   有了这些，用户那张表（22,006 / 29 / 0）可以从账本直接复算，而不必相信一段说明。
2. **保留降级前的判定**（`R/phenotype-evidence.R`）：identity demote 不再覆盖原值，
   新增 `<prefix>_call_pre_demote` 列与 `n_pre_demote_supported` / `n_pre_demote_probable`。
   报告因此可以同时呈现「证据说 29 个 supported」与「轴检验后全部 undetermined」。
3. **记录降级是数据集级的这件事**：降级语句旁的注释写明「一个数据集级统计量会移除所有肯定
   结论，使 malignant 不可达」，并要求报告同时呈现两态。

## 四、还没改、值得讨论的两点

1. **`transcriptomic_cut` 的默认是否该保留**。中位数默认让「签名+」永远约 50%，
   在报告里读起来像「一半细胞有签名」。可选方案：(a) 要求显式给 cut（拒绝静默默认）；
   (b) 保留默认但把状态改名（例如「签名高于池中位数」）并在报告里强制写出比例与来源。
   我倾向 (b) + 证据里的 `signature_cut_is_median` 标志——改默认会破坏既有工作流，
   而标志已经让读者无法误读。
2. **CNV 线该不该报「灵敏度未知」**。用户的数据说明：在 CNV 稀少的队列上，
   `non_malignant`（CNV− 显式否定）会把「看不到 CNV」写成「非恶性证据」——
   201 个 DNA 突变细胞里 90 个被判 non_malignant 就是这个失效模式。
   最小修法是在 malignant 报告里对 `cnv_reference_threshold` 与 CNV+ 比例加一条
   固定限制语句：「CNV 线为阴性不等于非恶性，除非该队列的 CNV 灵敏度已被验证」。

## 五、给回复的建议要点

1. 先承认两件事：签名线默认是**相对**判据（中位数），四态里的 `undetermined` 因此占大头；
   identity 降级是**数据集级**的，会让 `malignant` 在该队列上不可达——这不是他们的操作问题。
2. 指出 `hmm = FALSE` 使这条 CNV 线与包文档描述的默认路径不同，并说明我们已经把 `hmm`
   与 `cutoff` 记进证据（他们下次可以直接从账本复算）。
3. 对「201 个 DNA 突变细胞 0 个 CNV+」给出两层解释：阈值（P99+3MAD）确实保守；
   该队列可能 CNV 稀少。并同意他们的推论：这正说明**「两条线同时阳性」在 CNV 稀少的队列上
   几乎不可达**，四态的输出会整体偏向 undetermined / non_malignant。
4. 感谢 DNA 仲裁——它是这次最有价值的部分，因为它把「CNV− 被当作非恶性证据」这个失效模式
   量化了（90/201）。我们已经把「CNV− ≠ 非恶性」写成报告限制。
