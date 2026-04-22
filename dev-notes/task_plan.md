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

## Status
**In Progress** - 已完成 P0-P6 的状态层重构、统一读取/存在性接口建设，以及包级回归验证；当前正在推进 P7 的对外 API 命名统一。
