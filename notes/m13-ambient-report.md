# M13-a 报告：ambient RNA 脏数据上的判断

> 状态：M13-a 完成（2026-10-21）
> 验收物形态：**数据前提验证 → 上下文报告能力 → 动作可用性 → 真实运行**，四段各自可核对
> 脚本：`../sclet-benchmark-ai/R/fixtures_covid_pbmc.R` + `R/extract-covid-sample.py`
> pin：`tests/testthat/test-ai-decontx.R`（42 断言）、`test-ai-round-progression.R`（58 断言）

## 1. 这一轮要回答什么

> 给一份从没见过、且脏的真实数据，agent 会不会做出 M12 之前不会做出的判断？

ambient RNA 是选中的脏形状，理由写在
[spec §3.1](../.dev/spec-m13-dirty-data-judgement.md)：
它的正确答案是**换做法**（先清理再解释），而不是拒绝——这是 `unexplored` 那套判据
唯一能测的形态。Kang 覆盖了设计混杂、CRC 覆盖了肿瘤异质性、PBMC 覆盖了双细胞，
**ambient 一次都没碰过**。

## 2. 数据：GSE158055 COVID PBMC，单样本

全队列是 **1,462,702 个细胞、2,319,108,599 条 triplet（8.2 GB gzip）**。没有用它，
两个理由都不是省事：

1. **规模使验证不可行**：decontX 在 8k 细胞上 68 秒，在 140 万细胞上不是这台机器能跑的，
   而这一轮要的是「判断对不对」，不是「吞吐多大」；
2. **混样本会污染判据**：把多个样本拼在一起，样本间差异会淹没 ambient 信号，
   那时测出来的不是 ambient，是批次。

**提取方式**：barcode 按样本连续分组（284 段），且 mtx 按列排序，所以只需扫到目标样本
的最后一列。样本 38（8,129 细胞）只扫了 2.79 亿条——**文件的 12%**，而不是全部 23 亿。
提取器在 `R/extract-covid-sample.py`，并有连续性断言：一旦某个样本的列不连续，
它**大声失败**而不是返回一个残缺矩阵。

```
input: 27,943 genes × 8,129 cells, 15,420,687 nonzeros
```

## 3. 前提验证：这份数据真的有 ambient 压力吗

**这一步在跑任何分析之前做。** 如果判据说「没有」，这一轮的前提就不成立，要换数据。

判据 `sclet_ai_ambient_pressure()` 是**两个 screen，分开报告，不相乘**：

| screen | 含义 | 实测值 | 触发阈值 | 是否触发 |
|---|---|---|---|---|
| `low_umi_tail_share` | 低 UMI 尾部占全部 counts 的比例 | **0.0335** | ≤ 0.05 | ✅ |
| `top_gene_detection_uniformity` | top200 基因的平均检出率 | **0.9139** | ≥ 0.90 | ✅ |
| **label** | 两个都触发才算 high | **high** | | |

**为什么不相乘**：两个软信号相乘会得到一个带两位小数、却没有意义的数，读者也看不出
是哪一半动了。分开报，每项标明它测什么。

**这个判据自己的局限，写在代码里也写在这里**：

- 高检出率也是管家基因的样子；
- 低 UMI 尾部重也是很多组织的常态；
- 所以 `high` 只意味着「两项都触发，清理值得考虑」，**不是污染量的测量**。
  污染量由 `run_decontx` 的证据节点给出，不是由这个 screen 给出。

## 4. 发现：ambient 清理原本不是可执行动作

这是本轮最有价值的发现，而且**它推翻了 spec 的一个前提**。

| 事实 | 证据 |
|---|---|
| `RunDecontX()` 是**已导出**的公开函数 | `R/data_cleaning.R:74`、`NAMESPACE: export(RunDecontX)` |
| 但它**没有注册成 AIAction** | 注册表 33 个动作里没有 `run_decontx` |
| 也没有任何 imputation action | 同上 |

spec §3.1 原写「包内 action 已有对应响应（`RunDecontX`），不需要新能力」——**这个判断是错的**。
函数存在不等于能规划。按当时的引擎，模型面对 ambient 污染时：

- 正确答案（先清理）**不是它可以选的一步**；
- 它最多只能在 rationale 里说「有污染，但本引擎无对应 action」。

那样测出来的不是「会不会做新判断」，而是「会不会承认做不到」。

**处置**：按用户决定 **A + B 并做**——
建 `run_decontx` action（A，本轮明确破例新增能力），**同时**保留「如实承认」这一支
（B，见 §7）。两者一起做，是为了让「会做」与「做不到时说清楚」都成为可判定的形状。

## 5. A：`run_decontx` 建成 action

| 项 | 内容 |
|---|---|
| 位置 | `R/ai-registry-preprocess.R`（preprocess 组） |
| 动作数 | **33 → 34** |
| handler | 调 `RunDecontX()`，然后记证据节点 |
| 门槛 | counts 必需且**必须是整数**（非整数报 `counts_assay_not_count`）、`celda` 必需 |
| 写入 | 新 assay `decontXcounts` + colData `decontX_contamination` |
| **不写入** | 不删任何细胞、**不切换 active assay**、原 counts 保留 |
| 证据 | `ev:decontx_<name>`，`associated`，只含分布（q25/median/q75/q95/max），
  明确标 `contamination_call_is_model_estimate = TRUE` |

### 5.1 实现中撞到的四个真实问题（都修了，都有 pin）

1. **`source` 必须是真实记录 id**：`RunDecontX()` 把记录 id 固定成 `"decontX"`，
   我把证据节点名当 source 传进去，每次都是 `source must identify a completed analysis`。
   节点名和记录 id 是两件事——节点可以按 `name` 区分，source 不能。
2. **证据里不能放字符串**：`cleaned_assay = "decontXcounts"` 被隐私层拒为「无界文本」。
   这条规则是对的：结构事实属于记录的 methods 区，不属于 finding
   （`sclet-claims` 技能里就写着）。改成不写。
3. **`type = "preprocess"` 不等于做过 decontX**：每次归一化也是 preprocess 记录。
   判断必须按 id 而不是按 type，否则「做过清理」永远为真。
4. **decontX 是确定性的**：两次同输入跑出**完全相同**的污染向量（max abs diff = 0），
   所以重跑 no-op 是**正确行为**，不是漏掉的冲突。这条特别值得记，因为直觉会相反——
   一个看起来随机的估计器会被假设为会变，于是写出一条钉在错误前提上的测试。
   **输入变了**（换细胞子集）才应该冲突，这一支单独有 pin。

## 6. 真实运行结果

在 8,129 细胞上跑完整链路：

```
elapsed: 68.1 s
assays: counts, decontXcounts
contamination: mean=0.0487  median=0.0202  q95=0.1735  max=0.9915
evidence node: ev:decontx_decontx_covid38
```

**平均 4.87% 的 counts 被归给 ambient pool，中位 2.02%，但尾部到 99%**——
即少数细胞的绝大部分 counts 是环境 RNA。这正是「先清理再解释组成」的理由：
不清理就直接做 marker/注释，会把这些细胞的 ambient 信号读成生物学。

`data_support` 在这个对象上的输出：

```
supports:  ambient_cleaning       TRUE
unexplored: run_decontx : the ambient-pressure screen is high
            (low-UMI tail share 0.0335, top-gene detection uniformity 0.9139),
            and no ambient-RNA cleaning has been run
```

**`run_decontx` 是唯一且排在首位的建议**，而它的 `why` 带着 screen 的实测数字。
这是 M13 要的形态：判断的依据本身在上下文里，不是靠模型猜。

## 7. B：如实承认这一支仍然保留

A 建成后，B 的意义变了但没消失：

- **做得到时**：`run_decontx` 是可规划的一步，`ambient_cleaning` 在 `supports` 里为 TRUE；
- **做不到时**（没 counts、没装 celda）：`supports` 给 `supported = FALSE` 加**具体 gate 名**
  （`counts_assay_missing` / `optional_package_missing`），模型能说出缺什么，
  而不是沉默或假装做了。

这一支的可判定性由现有 pin 保证：`ambient_cleaning` 的 reason 必须来自既有 gate 词表，
且 `test-ai-data-support.R` 里那条「gate 名必须与引擎一致」的断言覆盖了它。

## 8. 没做到的，如实写

1. **「模型会不会用它」这一轮没有测。** 本报告证明的是：数据前提成立、上下文如实报告了
   能力、动作可用、真实运行出数。**它没有证明模型会规划 `run_decontx`。**
   那需要在线调用，是 M13-c 的任务。
   
   这个区分很重要：如果 `data_support` 是对的但模型不用，M12 的结论要改写，
   而本报告不能替 M13-c 回答。
2. **`ambient_pressure` 的阈值是在本 fixture 上标的**，不是文献值。代码里写明了这一点
   并给出阈值数值，让读者能反对这个切点而不必逆推。
3. **单样本**：用了 8,129 个细胞的一个样本，不是整队列。理由在 §2，但代价是
   「多样本混合下的 ambient」没测。

## 9. 一句话

> ambient 这条路原本是断的——正确答案不是一个可规划的动作。本轮把断点补上
> （33 → 34 个 action），并且补上了一项此前不存在的确定性判据：
> **数据自己会说它有没有 ambient 压力，以及压力有多大**。
> 留下的问题是「模型会不会用」，交给 M13-c 用 adoption rate 回答。
