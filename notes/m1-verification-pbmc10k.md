# M1 复验记录（PBMC 10k 报告）

> 独立进程复验：`readRDS("m1-report.rds")` 后只读账本重算，不采信生成脚本的任何中间值。时间 2026-10-20 14:20:53。
> 账本指纹 `sclet-ai-593d0429`；对象 33538 × 11769。

| # | 断言 | 结果 | 细节 |
|---|---|---|---|
| V1 | 记入账本的报告存在且含 4 条结论 | PASS | findings: 4 |
| V2 | 所有引用的证据 id 都能解析 | PASS | refs: 3 |
| V3 | 每条记录的等级 = 引擎重算的等级 | PASS | 4/4 |
| V4 | 没有结论超过其引用的证据上限 | PASS | 0 overclaims |
| V5 | 报告里的 11 个关键数字都能在账本里重算并在正文出现 | PASS | matched: 11 |
| V6 | 报告引用的 ceiling 结论可重算（markers + GO = 1 条证据线） | PASS | ceiling: associated groups: 1 |
| V7 | 聚类 seed 与 resolution 记录在账本里 | PASS | seed: 1 resolution: 0.8 |
| V8 | 证据节点里没有富集条目名/基因名泄漏 | PASS | numeric/aggregate only |
| V9 | 报告里写的账本指纹 = 记录时存下的 context_fingerprint | PASS | report: sclet-ai-4d3dfbea recorded: sclet-ai-4d3dfbea |
| V10 | markers + GO + KEGG 只算一条独立证据线 | PASS | groups: 1 ceiling: associated |
| V11 | RenderAIReport() 从账本渲染并通过自检 | PASS | self-check passed |

总计：11/11 通过。
