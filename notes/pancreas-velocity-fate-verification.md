# M11 复验记录（真实数据 速度-轨迹-命运 报告）

> 独立进程只读账本重算。时间 2026-10-20 23:30:58。
> 账本指纹 `sclet-ai-54e8042d`；对象 27998 × 3696。

| # | 断言 | 结果 | 细节 |
|---|---|---|---|
| V1 | 报告存在且含 4 条结论 | PASS | findings: 4 |
| V2 | 引用的证据 id 可解析 | PASS | refs: 3 |
| V3 | 每条记录等级 = 引擎重算等级 | PASS | 4/4 |
| V4 | 对象确有 spliced/unspliced/counts 三套 assay（速度的真实性前提） | PASS | counts,spliced,unspliced,logcounts |
| V5 | 去向分布 = 从吸收概率列重算的 argmax 分布 | PASS | recomputed: 0.1071/0.599/0.2938 |
| V6 | terminal_state_fractions 只统计终态细胞（分母写进节点） | PASS | terminal cells: 90 of 3696 |
| V7 | lineage 数 = trajectory 记录 | PASS | lineages: 5 |
| V8 | 报告区分了终态占比与去向分布 | PASS | both fields named |
| V9 | 报告没有把匿名终态编码写成 alpha/beta/delta | PASS | anonymized codes stated |
| V10 | RenderAIReport() 重渲染并通过自检 | PASS | self-check passed |

总计：10/10 通过。
