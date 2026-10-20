# M5 复验记录（PBMC 10k 注释报告）

> 独立进程复验：只读账本重算，不采信生成脚本的中间值。时间 2026-10-20 15:32:39。
> 账本指纹 `sclet-ai-52685d27`；对象 33538 × 11769。

| # | 断言 | 结果 | 细节 |
|---|---|---|---|
| V1 | 报告存在且含 4 条结论 | PASS | findings: 4 |
| V2 | 所有引用的证据 id 可解析 | PASS | refs: 2 |
| V3 | 每条记录等级 = 引擎重算等级 | PASS | 4/4 |
| V4 | 报告里的主导群数/混合群数/分数/参考名都能重算并出现 | PASS | matched: 4 |
| V5 | 每个主导标签都出现在报告的组成行里 | PASS | labels: 7 |
| V6 | 参考来自已记录的 user_decision | PASS | decisions: 1 |
| V7 | markers + annotation = 2 条独立线，上限 associated | PASS | groups: 2 ceiling: associated |
| V8 | RenderAIReport() 重渲染并通过自检 | PASS | self-check passed |
| V9 | 混合群没有被写成确定的细胞类型 | PASS | mixed clusters: 9 |

总计：9/9 通过。
