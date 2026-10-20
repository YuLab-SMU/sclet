# M8 复验记录（PBMC 10k 细胞通讯报告）

> 独立进程只读账本重算；互作数经 get_cci() 表重算。时间 2026-10-20 17:59:04。
> 账本指纹 `sclet-ai-50126a6b`；对象 33538 × 11769。

| # | 断言 | 结果 | 细节 |
|---|---|---|---|
| V1 | 报告存在且含 3 条结论 | PASS | findings: 3 |
| V2 | 引用的证据 id 可解析 | PASS | refs: 2 |
| V3 | 每条记录等级 = 引擎重算等级 | PASS | 3/3 |
| V4 | 证据节点的互作数/发送群/接收群 = get_cci() 表重算 | PASS | evidence: 733/12/13 | recomputed: 733/12/13 |
| V5 | 引擎来自已记录的 user_decision（cci_method = CellChat） | PASS | decisions: 2 |
| V6 | 互作数/发送群/接收群/概率中位数都能重算并出现 | PASS | matched: 4 |
| V7 | 报告写明 cluster 12（双细胞群）参与了互作推断 | PASS | cluster 12 doublet fraction: 0.983 |
| V8 | 证据节点不含配体/受体/通路名称 | PASS | counts and quantiles only |
| V9 | RenderAIReport() 重渲染并通过自检 | PASS | self-check passed |

总计：9/9 通过。
