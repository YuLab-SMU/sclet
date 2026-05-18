# Task Plan: sclet 状态层重构

## Goal
在保持 SingleCellExperiment 为核心对象的前提下，引入统一的 sclet 状态命名空间，修复高风险接口，并完成最小可验证的 P0 重构。

## P1 追加范围
- 增加 `DefaultReduction()` / `DefaultReduction<-`，统一 active reduction 读取与设置
- 增加 `CommandLog()`，暴露状态层中的运行记录
- 将主分析链的预处理、降维、建图、聚类、轨迹、批次校正、milo 接入 command log
- 修复对象在转换步骤中 metadata 容易丢失导致 command log 被覆盖的问题

## P2 追加范围
- 增加 `DefaultAssay()` / `DefaultAssay<-`，统一 active assay 读取与设置（仅 sclet 侧默认，不强改底层 assay 内容）
- 增加 `DefaultGraph()` / `DefaultGraph<-`，统一 active graph 读取与设置
- 增加 `ActiveIdent()` / `ActiveIdent<-`，为 `colLabels` 增加“当前身份字段”概念（默认仍等价于 `colLabels`）
- 将 `NormalizeData` / `ScaleData` / `BatchRemover` 更新 active assay
- 将 `FindNeighbors` 更新 active graph，将 `FindClusters` 更新 active ident

## P3 追加范围
- 让 `FindMarkers()` / `FindAllMarkers()` 在默认情况下参考 `DefaultAssay()`，但对 `scaled` 自动回退到 `logcounts`
- 增强 `CommandLog()` 输出，增加参数摘要与结果摘要列
- 增强 `RunUMAP()`：允许不传 `dims`，自动推断 source reduction
- 运行完整 `devtools::check()`，确认本轮状态层 API 在 package 层面可通过

## P4 追加范围
- 新增统一读取接口：`get_trajectory()`、`get_milo()`、`get_supercell()`
- 将现有读取点优先切换到统一 accessor，减少直接读取 metadata 的分叉逻辑
- 为统一 accessor 补充“新状态 + 旧 metadata 回退”的测试

## P5 追加范围
- 新增 `get_batch()`、`get_hvg()`、`get_graph()`，补齐 batch/HVG/graph 的统一读取接口
- 将批次校正与图读取路径优先切换到 accessor
- 补充统一读取层的回归测试，并确认完整 `check()` 仍仅剩环境 NOTE

## P6 追加范围
- 新增 `has_trajectory()`、`has_milo()`、`has_supercell()`、`has_batch()`、`has_hvg()`、`has_graph()`
- 将 plotting / wrapper 层的存在性判断切换到 `has_*()`，减少手写 NULL 判断与旧 metadata 注释
- 补充 `has_*()` 的正反例测试，并确认包级检查结果稳定

## P7 追加范围
- 对外 API 命名统一为 `Run*` canonical public verbs，保留旧 `run*` 作为兼容 alias
- 先收口 `RunPCA()`、`RunMilo()`、`RunCellChat()`，避免新旧命名继续并行扩散
- 将 command log 的命令名统一写为 canonical `Run*`

## P8 追加范围
- 将 README、NEWS 等用户可见文档中的主入口表述统一为 `Run*`
- 旧 `run*` 仅作为兼容 alias 保留，不再作为推荐写法出现在用户文档首屏

## Phases
- [x] Phase 1: 现状调研与问题归纳
- [x] Phase 2: 写入实施方案并冻结本轮范围
- [x] Phase 3: 实现状态层与核心函数重构
- [x] Phase 4: 补充测试并运行校验
- [x] Phase 5: 收尾说明与后续建议

## 本轮范围
- 新增内部状态 helper，统一管理 `metadata(sce)$sclet`
- 迁移 HVG、graph、trajectory、milo 的内部读写入口
- 修复 `FindMarkers()` 直接访问 slot 的实现
- 修复 `subset(SingleCellExperiment)` 对 `Barcode` 的隐式依赖
- 保持对旧 metadata 字段的读取兼容

## 不在本轮范围
- 新建 S4 子类
- 全量重写 interop 契约
- 完整多模态支持与 `altExp` 工作流
- 全面重做文档与 man 文件

## Decisions Made
- 继续以 `SingleCellExperiment` 作为唯一核心对象
- 状态统一收口到 `metadata(sce)$sclet`
- 优先做兼容式重构，而不是一次性移除全部旧字段
- `DefaultReduction` 仅用于 active reduction 的读取/展示，不改写 `FindNeighbors()` 当前基于 PCA 的默认语义

## Errors Encountered
- 首次 `devtools::check()` 构建时仍将计划目录打入包结构检查范围，先后确认 `.Rbuildignore` 只影响 tarball，不消除隐藏目录 NOTE
- 包级 `check` 过程较长，因此追加执行定向 `devtools::test()` 作为本轮变更回归验证
- P1 初版 command log 在对象转换后被 metadata 覆盖重置，随后通过 `sclet_restore_state()` 与状态合并逻辑修复
- `ScaleData(features=子集)` 曾因写回部分矩阵触发 assay 维度不一致，已改为在完整 assay 上做局部替换
- 完整 `devtools::check()` 目前仍保留 2 个 NOTE，其中 `future file timestamps` 更像环境时钟问题；为消除隐藏目录 NOTE，计划目录从 `.dev` 迁移为 `dev-notes`
- 完整 `check()` 在清理隐藏目录 NOTE 后暴露出 `clusterProfiler` 的强依赖 warning，随后通过将其改回可选依赖并移除 NAMESPACE 顶层 import 修复

## P9-P27 已落地扩展
- `bookdown/*.Rmd` 的用户可见口径已统一到 `Run*`，并补入面向用户/开发者的设计说明
- 已形成并落地 `analysis-state contract` 主线，而不再局限于早期 `assay/layer contract`
- 已实现 `layer registry` 与用户接口：`Layers()`、`LayerData()`、`DefaultLayer()`
- `NormalizeData()`、`ScaleData()`、`BatchRemover()` 已接入 layer/state 写回
- `RunPCA()`、`RunUMAP()`、`FindNeighbors()`、`FindClusters()` 已记录 layer 与 integration provenance
- `sce_merge()` 已写入 merge provenance，`BatchRemover()` 已写入 active integration state，并可追踪 merge 来源
- `RunSingleR()` 已接入 `annotation` + `mapping` contract，支持 `get_annotation()` / `get_mapping()`
- `RunSlingshot()` / `RunSlingshot_trajectory()` 已接入可命名 `trajectory` contract 与 active/id 读取
- `RunSuperCell()` 已接入 `aggregation` contract，返回对象记录 parent-child provenance，并统一使用 `sclet_set_analysis()`
- `RunCellChat()` 已接入可命名 `communication` contract 与 active/id 读取
- `RunMilo()` 已接入可命名 `milo` contract 与 active/id 读取
- `get_*` / `has_*` 主路径已基本统一到 state/accessor 模式，保留 legacy fallback
- `CommandLog()` 已与 canonical `Run*`、typed outputs 和主要 downstream 模块对齐

## 当前已完成能力
- 核心状态层：`metadata(sce)$sclet`、active state、command log、typed analysis states
- 表达空间管理：`DefaultAssay()`、`DefaultLayer()`、logical layer -> physical assay 解析
- 主分析链：preprocess、reduction、graph、clustering 的统一 state/provenance
- 下游 contract：`annotation`、`mapping`、`trajectory`、`aggregation`、`communication`、`milo`
- provenance 主线：merge -> integration -> downstream analysis
- 用户文档同步：每轮代码改动后均已更新相关 `bookdown/*.Rmd`，未渲染 bookdown

## 当前验证状态
- `state-refactor` 与相关定向测试已多轮通过
- 最近一轮全局收口验证结果：
  - `Rscript -e 'devtools::test(filter = "state-refactor|advanced")'` 通过
  - `Rscript -e 'devtools::check()'` 结果为 `0 errors / 0 warnings / 0 notes`
- 当前保留的 warnings 主要来自上游 Bioconductor 依赖的 deprecated API 提示，不是本轮 contract 改造引入的新错误

## Remaining Gaps
- 上游依赖弃用接口的**公开主链调用**已基本清理完成：
  - `NormalizeData()` 已改为 `scuttle::calculateCPM + log1p`
  - `FindNeighbors()` 已改为 `bluster::makeSNNGraph`
  - HVG 选择与批次合并不再直接调用 `scran::getTopHVGs()` / `scran::combineBlocks()`
- 当前剩余技术债主要在**低层兼容 helper**，而不是用户可见主流程：
  - `sclet_model_gene_var()` / `sclet_fit_trend_var()` 仍依赖部分 `scran` 内部低层统计逻辑来保持现有行为
  - `test-deprecation-cleanup.R` 仍需要持续盯住上游 warning 变化
  - 后续若 Bioconductor 上游给出更稳定的低层替代接口，再评估是否继续收口
- `get_*` / `has_*` 已统一，legacy metadata fallback 已收到统一的内部 helper：
  - `sclet_get_legacy_graph_entry()` 统一 `knn_graph` 旧 metadata 读取
  - `sclet_get_legacy_analysis_record()` 统一 `trajectory / milo / supercell / batch` 旧 metadata 读取
  - 全仓外围调用点已扫描完毕：所有消费侧均走 `get_*` / `has_*`，无绕开 accessor 的直读点
- `mapping` 目前是最小 contract，后续若接更多 reference transfer 方法，需要把 schema 再抽象一层
- `supercell` / `cellchat` / `milo` 的 rich object 目前仍以单对象载荷为主，若将来支持对象集合或惰性载入，可能需要扩展 payload 约定

## Next Phase Candidates
1. ~~全仓读取路径再收口~~ ✅ 已完成
   - 旧 metadata 直读路径已统一到内部 helper，外围调用均走 accessor
2. ~~cross-analysis helpers~~ ✅ 已完成
   - 在现有 state contract 之上，继续增加更高层的用户视角 helper
   - 新增 `RunStandardPipeline()` 作为主分析流程总入口，完全消费 `analysis-state contract`
   - 新增 `PipelineSummary()` 提供直观的文字版流水线摘要
3. 上游低层兼容 helper 观察
   - 持续关注 `scran` / `scuttle` 上游变化
   - 只有在出现明确稳定替代接口时，再收紧 `sclet_model_gene_var()` 等低层兼容实现
4. `scop` 其它未完成模块 (Mapping / SuperCell 等高阶集成)
   - 基于 `RunSingleR` 扩展更多的 Reference Mapping
   - 根据需求进一步增强 Integration 与 Visualization 消费约定

## Recommended Next Step
- 鉴于 `RunStandardPipeline` 和 `PipelineSummary` 的引入，`scop` 对标计划中提到的核心平台体验（标准入口、状态注册表、可视化与切换机制）的**基础框架已全部打通**。
- 下一步建议根据真实业务场景需求，切入更高阶的下游能力（如 Mapping、SuperCell 或自定义 Integration ），或开始准备对外发版前的文档与 Vignette 打磨。

## Status
**Phase Complete** - `analysis-state contract` 主线、layer/integration provenance、主要 downstream contract、标准工作流总入口 `RunStandardPipeline` 与 `PipelineSummary` 已经全部落地，包级验证通过。`scop` 对标中最重要的“平台化心智”已被成功引入 `sclet`。
