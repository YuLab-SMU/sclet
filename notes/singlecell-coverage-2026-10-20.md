# 单细胞覆盖盘点（截至 2026-10-20）

> 回答一个问题：**“单细胞这部分是不是已经看过一遍了？”**
> 答案取决于“看过”的标准，所以这里把两个标准分开列。数据来源：注册表实测
> （33 个动作）、`../sclet-benchmark-ai/docs/action-coverage.md`（动作级记录）、
> `notes/m1-report-*` / `notes/m2-report-*` / `notes/pbmc10k-annotation-*`（闭环记录）。

## 两个标准

**A. 动作级**：某个动作被模型在真实数据上规划并执行过至少一次，有 benchmark 记录。
覆盖表 2026-10-19 重数后，除 `run_milo`、`inspect_ledger` 外都有记录；
spatial/multimodal 六个是 typed refusal（数据不具备，并非未尝试）。

**B. 闭环级**（M1/M2 的标准）：真实数据 → 证据节点 → 引擎评级 → `RenderAIReport()`
报告 → 独立进程只读账本复验。这一档才算“看过一遍”。

## 当前状态（33 个动作）

| 域 | 动作 | 动作级 | 闭环 |
|---|---|---|---|
| 主干 pipeline | normalize / HVG / scale / PCA / neighbors / clusters | run | ✅ M1（PBMC 10k，43 min 内跑完，含报告） |
| 标记基因 | `run_de_test`（全群 + 条件配对） | run | ✅ M1 / M2（pseudobulk limma，16 样本） |
| 富集 | `run_enrichment`（GO + KEGG） | run | ✅ M1（6,430 / 645 条） |
| 批次整合 | `run_integration` / `compare_integration_routes` / `confirm_design` | run | ✅ M2（Kang，3 路线 + scVI 子集） |
| 细胞类型注释 | `run_annotation` | run | ✅ M5（MonacoImmuneData，13 群 → 7 类） |
| doublet / rare | `run_doublet_detection` / `run_rare_cell_detection` | run | ❌ |
| 轨迹 / 速度 / fate | `run_trajectory` / `run_velocity` / `run_fate_analysis` / `confirm_fate_terminal_states` / `summarize_fate_drivers` | run / gated | ❌ |
| 细胞通讯 | `run_cci` | run（gated） | ❌ |
| 差异丰度 | `run_milo` | run | ✅ M6（Kang，25 units，3,414 邻域） |
| 拷贝数 / 恶性 | `run_cnv_inference` / `call_malignant_cells` | run（GSE144735） | ❌ |
| 检查类 | `check_qc` / `inspect_status` | run | ✅（M1 流水线内） |
| 账本只读 | `inspect_ledger` | **none** | ❌ |
| 空间 / 多模态 | 6 个只读门 | gated（拒） | 不属纯单细胞范围 |

## 仍然缺的（纯单细胞）

1. **doublet / rare**：`smallest_group_n = 39` 这类数字只有配上 doublet 证据才有意义。
2. **trajectory → velocity → fate**（一条链，别拆三次）。
3. **`run_milo`**：唯一连动作级都没跑过的单细胞动作；Kang 的 donor×condition 正好能做。
4. **`run_cci`**：需要一次带用户决策（引擎选择）的闭环。
5. **CNV / malignant**：最贵（实测 66 min、22.6 GB），有 GSE144735 记录，闭环可最后做。
6. **`inspect_ledger`**：便宜，顺手。

## 数据形状的缺口

已闭环的两次都是“顺路”数据：PBMC 10k 是单样本，Kang 是规整的双 batch + 8 配对 donor。
真实单细胞的难点——多样本多条件、肿瘤异质性、ambient RNA、超稀疏、批次与条件深度混杂——
只碰到了最后一项（Kang 的 stim 只在 batch2），报告里已标为警告。
所以更准确的说法是：**主干成立，两种最常见的数据形状能走通闭环**，而不是“单细胞看完了”。

## 剩余工作量估计

每条闭环按 M1（半天）/ M2（半天，含途中挖出的 claim 级度量缺陷）/ M5（半天，含澄清通道）
的节奏，5–6 条大约 **5–7 天**。别按纯叠加计划：M2 那种“顺手发现度量测错对象”的比例是常态。
