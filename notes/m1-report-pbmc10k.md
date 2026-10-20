# PBMC 10k 端到端报告（M1 验收物）

> 数据：10x PBMC 10k，33538 基因 × 11769 细胞（只有 counts）。
> 路径：pipeline → markers → GO → KEGG，四轮分轮计划，每轮重新规划并重新校验。
> 账本指纹（记录本报告之前）：`sclet-ai-4d3dfbea`；记录后 `sclet-ai-281bee62`。
> 报告记录 id：`ai_report_pbmc10k_m1`；sclet 1.0.2 @ `9400e1d`；生成时间 2026-10-20 13:39。

## 方法与参数（描述性，不由引擎评级）

| 环节 | 实现 | 关键参数 |
|---|---|---|
| 标准化 | scuttle::calculateCPM + log1p | scale.factor = 10000 |
| 高变基因 | scran（引擎默认方法） | nfeatures = 2000 |
| PCA / 近邻 | scater::runPCA / bluster::makeSNNGraph | ncomponents = 50, dims = 1:30 |
| 聚类 | igraph::cluster_louvain | resolution = 0.8, seed = 1 |
| 标记基因 | FindAllMarkers | min.pct = 0.1, logFC >= 0.25 |
| 富集（GO） | clusterProfiler::compareCluster(enrichGO) | db = GO, ont = BP, keyType = SYMBOL |
| 富集（KEGG） | clusterProfiler::compareCluster(enrichKEGG) | db = KEGG, keyType = ENTREZID（表由 bitr 转换） |

## 结论（等级由引擎按证据判定）

| 结论 | 证据引用 | 记录等级 | 引擎审计 |
|---|---|---|---|
| 在 resolution 0.8、seed 1 下得到 13 个 louvain 群；对全部群做标记基因检验，共 12588 条基因-群记录（min.pct 0.1，logFC >= 0.25，最小群 39 个细胞）。 | `ev:findallmarkers` | associated | engine: ok |
| marker 基因在 GO BP 上得到 6430 条富集条目，覆盖 13 个群。 | `ev:enrichment_enrich_go` | associated | engine: ok |
| 同一批 marker 转成 ENTREZID 后在 KEGG 上得到 645 条富集条目，覆盖 13 个群。 | `ev:enrichment_enrich_kegg` | associated | engine: ok |
| 这些群可能对应不同的 PBMC 细胞类型（例如含一个约 39 个细胞的小群），但本次没有做基于参考的注释。 | — | hypothesis | engine: ok |

### 为什么最高只到 `associated`

`check_claim_ceiling()` 对这组证据给出：`associated`（all supporting evidence resolves to 1 independent line(s) of support; a claim cannot exceed 'associated'）。

GO 与 KEGG 是同一份 marker 表的两种注释，`ev:enrichment_*` 的 `parents` 指向 `ev:findallmarkers`，
所以它们和 marker 检验算**一条**证据线；换一个数据库不构成第二次独立测量。要把等级推到
`consistent_with`，需要一条真正独立的证据（例如基于参考的注释，或另一个样本/模态）。

## 明确没有说的话

- 没有说任何 `cluster_N` 是某类细胞（未做参考注释，相关结论为 `hypothesis`）。
- 没有说任何条件/批次效应（单样本，无 design）。
- 没有引用富集条目名（证据节点按设计只存计数与分位数；名称在 `get_enrichment()` 的状态记录里）。

## 复验方式

1. `readRDS(m1-report.rds)`，只读账本重算每个数字与等级，不读本文。
2. 断言：每个引用的证据 id 可解析；每条记录等级 = 引擎重算等级；每个数字与账本记录相等。
复验记录见 `notes/m1-verification-pbmc10k.md`。

## 下一步（建议，未执行）

- 若需要细胞类型，用用户指定的 celldex 参考跑 run_annotation。
- 若要评价条件效应，先确认 design 角色，再用 run_de_test(group=, covariates=)。

> 警告：只有一个样本，没有条件对比，也没有 donor/batch 结构，因此不能做条件或批次层结论。
> 警告：没有做基于参考的细胞类型注释（需要用户指定 celldex 参考）。
