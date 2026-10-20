# M7 复验记录（PBMC 10k 双细胞 / 稀有群报告）

> 独立进程只读账本重算。时间 2026-10-20 17:15:48。
> 账本指纹 `sclet-ai-156a1131`；对象 33538 × 11769。

| # | 断言 | 结果 | 细节 |
|---|---|---|---|
| V1 | 报告存在且含 4 条结论 | PASS | findings: 4 |
| V2 | 引用的证据 id 可解析 | PASS | refs: 3 |
| V3 | 每条记录等级 = 引擎重算等级 | PASS | 4/4 |
| V4 | 证据节点的细胞数/双细胞数/比例 = 从 colData 重算 | PASS | evidence: 863/11769 (0.073328) | recomputed: 863/11769 (0.073328) |
| V5 | 全局比例/最高群比例/cluster 9 比例/稀有群数/极大值数都能重算并出现 | PASS | matched: 4 | worst cluster: 12 |
| V6 | 零个稀有群时不注册证据节点，报告把它写成 hypothesis | PASS | rare clusters: 0 | rare nodes: 0 |
| V7 | 没有细胞被移除（对象细胞数 = 注释时的细胞数） | PASS | cells: 11769 | classes: doublet,singlet |
| V8 | RenderAIReport() 重渲染并通过自检 | PASS | self-check passed |

总计：8/8 通过。
