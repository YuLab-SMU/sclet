# M6 复验记录（Kang milo 差异丰度报告）

> 独立进程只读账本重算；置换对照读它自己的对象。时间 2026-10-20 16:21:50。
> 账本指纹 `sclet-ai-1585ff59`；对象 32738 × 43095。

| # | 断言 | 结果 | 细节 |
|---|---|---|---|
| V1 | 报告存在且含 4 条结论 | PASS | findings: 4 |
| V2 | 引用的证据 id 可解析 | PASS | refs: 2 |
| V3 | 每条记录等级 = 引擎重算等级 | PASS | 4/4 |
| V4 | 证据节点的邻域数/显著数 = 从 da_results 重算的值 | PASS | evidence: 1821 / 3414 | recomputed: 1821 / 3414 |
| V5 | 报告里的邻域数/显著数/中位 p/置换结果/DE 行数/对调 donor 数都能重算并出现 | PASS | matched: 7 |
| V6 | 置换标签的对照显著数为 0 | PASS | null: 0 / 3414 |
| V7 | 记录里能读到未配对的设计公式，报告据此写下限制 | PASS | formula: ~0 + condition + batch |
| V8 | RenderAIReport() 重渲染并通过自检 | PASS | self-check passed |

总计：8/8 通过。
