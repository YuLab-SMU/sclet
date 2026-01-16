# sclet 开发计划（里程碑驱动）

本计划的目标是把 `sclet` 提升为“完全可用、功能齐全、可维护”的单细胞分析工具包：以 `SingleCellExperiment` 为核心数据结构，提供稳定的 Seurat 风格入口函数、完善的多样本分析能力，并形成可复现的分析闭环。

## 设计原则

- **以 SCE 生态为主**：优先使用 `scater/scran/scuttle/SingleCellExperiment` 等 Bioconductor 生态能力。
- **可视化优先 ggsc**：所有通用可视化接口优先以 `ggsc` 实现；如果 `ggsc` 暂无对应能力，则仅记录需求并纳入 `ggsc` 开发计划，当前阶段不在 `sclet` 内实现。
- **依赖分层**：核心流程不依赖 Suggests；可选模块（如 CellChat/Slingshot/SuperCell/ComplexHeatmap）必须显式检测安装状态并给出清晰错误。
- **对象约定清晰**：统一细胞身份（ident）、批次（batch）、样本（sample）、细胞 ID（cell_id）等字段的读写约定，保证全链路一致。
- **可复现与可审计**：关键步骤的参数、随机种子、版本信息、输入摘要写入 `metadata()`，支持结果追溯。

## Milestone 工作流（每次完成一个里程碑都执行）

1. 提升版本号（DESCRIPTION），并在 NEWS.md 增加同版本段落。
2. 同步更新文档：man、README、bookdown。
3. 补齐 tests：覆盖新增/改动的 build 与渲染入口。
4. 质量门禁：
   - `Rscript -e 'devtools::document()'`
   - `Rscript -e 'devtools::test()'`
   - `Rscript -e 'devtools::check()'`
5. git 提交：把版本号、NEWS、tests、文档一起 commit。
6. 进入下一个里程碑：从 dev-plan.md 的下一段开始实现，完成后重复流程。

## Milestone 0.1.0：核心流程可用与一致性修复 (Completed)

### 目标

- 打通从导入到聚类与 marker 的最短闭环，并保证对象字段约定一致。

### 交付物

- 统一 ident 读写规则：聚类结果与下游 marker/可视化对齐（Idents/FindClusters/FindAllMarkers 全链路一致）。
- 修复阻断性 bug（例如 subset、基因注释汇总等）并补齐最小回归测试。
- 所有 Suggests 模块增加依赖检测与错误信息（CellChat/Slingshot/SuperCell/batchelor/ComplexHeatmap 等）。
- 核心预处理尽量“去 Seurat 内部实现依赖”：默认路径优先 `scuttle::logNormCounts()` 等。

### 验收标准

- 在示例数据上可复现完成：Read10X → QC → Normalize → HVG → PCA/UMAP → Neighbors → Clusters → Markers。
- `devtools::check()` 通过。

### 可视化策略

- 仅整理可视化入口与对象字段约定，不新增大量绘图函数。
- 涉及降维散点的展示优先使用 `ggsc::sc_dim()`；其他缺口统一写入需求文档（见 `ggsc-visualization-requirements.md` 与 bookdown 对应章节）。

## Milestone 0.2.0：QC 与可视化“日常可用”（ggsc 优先）(Completed)

### 目标

- 用户无需写额外 glue 代码即可完成常规 QC、探索性可视化与基础标注输出。

### 交付物

- QC 套件：标准 QC 指标展示与一键过滤（SCE 版），并可与后续聚类/注释衔接。
- 可视化入口统一：DimPlot/FeaturePlot/DotPlot/Heatmap/RidgePlot 等同名 API 的对齐方案。
- 注释入口（非可视化）：封装参考映射/marker-based 注释的最小闭环（输出落入 `colData`）。

### 验收标准

- 文档（bookdown）给出从 QC 到注释的可运行示例。
- 所有新增可视化如 `ggsc` 已支持则直接调用；不支持的仅在需求文档中列出，不在本里程碑实现。

### 说明

- `ggsc` 缺口统一记录，当前阶段不在 `sclet` 内另起实现，避免重复造轮子。

## Milestone 0.3.0：多样本与实验设计（可信差异分析）(Completed)

### 目标

- 从“单样本探索”升级为“多样本、可重复、统计可信”的分析能力。

### 交付物

- 增加 pseudobulk 差异分析工作流（按样本聚合后进行 DE），替代仅逐细胞 wilcox 的默认叙事。
- 批次整合策略矩阵：统一入口与输出规范（fastMNN 等为起点），并记录参数到 `metadata()`。
- 结果导出：标准化输出 marker/DE 表格、注释表、降维坐标与关键图表。

### 验收标准

- 多样本示例可以复现实验设计与差异分析。
- `devtools::check()` 通过。

## Milestone 1.0.0：工程化、互操作与性能 (Completed)

### 目标

- API 稳定、测试充分、与 Seurat/AnnData 互操作顺滑，可用于长期维护与生产分析。

### 交付物

- 互操作：SCE ↔ Seurat ↔ AnnData 的导入导出与字段映射规范。
- 性能：大对象（DelayedArray/磁盘后端）友好、关键步骤并行化策略。
- 质量体系：覆盖主流程的回归测试与持续集成策略。

### 验收标准

- 核心 API 变更有清晰 deprecation 策略。

## Milestone 1.1.0：高级分析模块与交互式探索 (Planned)

### 目标

- 扩展 `sclet` 的分析边界，涵盖轨迹推断、富集分析及交互式数据探索，提供一站式解决方案。

### 交付物

1.  **轨迹推断封装**：
    - 实现 `RunSlingshot`：封装 `slingshot` 流程，支持从降维结果自动推断谱系。
    - 提供简单的轨迹可视化入口（依赖 `ggsc` 或原生绘图）。
2.  **富集分析集成**：
    - 实现 `RunEnrichment`：无缝衔接 `FindAllMarkers` 结果与 `clusterProfiler`，一键进行 GO/KEGG 分析。
3.  **交互式探索**：
    - 引入 `iSEE` 依赖，实现 `RunExplorer()`，启动交互式 Shiny 应用，方便无代码背景用户探索数据。
4.  **Metacell 支持**：
    - 封装 `SuperCell` 流程 (`RunSuperCell`)，支持大规模数据的 Metacell 构建与下游分析。

### 验收标准

- 新增模块在 bookdown 中有完整演示。
- 交互式应用能正常启动并展示核心图表。
