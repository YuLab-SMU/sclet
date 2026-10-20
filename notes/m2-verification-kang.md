# M2 复验记录（Kang IFN-beta 报告）

> 独立进程复验：`readRDS("m2-report.rds")` 后只读账本重算，不采信生成脚本的中间值。时间 2026-10-20 10:40:00。
> 账本指纹 `sclet-ai-07216caf`；对象 32738 × 43095。

| # | 断言 | 结果 | 细节 |
|---|---|---|---|
| V1 | 记入账本的报告存在且含 4 条结论 | PASS | findings: 4 |
| V2 | 所有引用的证据 id 都能解析 | PASS | refs: 2 |
| V3 | 每条记录的等级 = 引擎重算的等级 | PASS | 4/4 |
| V4 | 没有结论超过其引用的证据上限 | PASS | 0 overclaims |
| V5 | 三条路线的 batch_mixing 与重算值一致（各在自身嵌入上） | PASS | routes: 3 |
| V6 | 报告里的关键数字（维度/基线 mixing/DE 行数/显著数/样本数）都能重算并出现 | PASS | matched: 6 |
| V7 | 配对设计来自 DE 记录而非文本（~ condition + donor, 16 样本） | PASS | ~ condition + donor | 16 samples |
| V8 | 比较证据节点不含路线名（只有 route_N 编码） | PASS | codes: route_1,route_2,route_3 |
| V9 | 报告里写的账本指纹 = 记录时存下的 context_fingerprint | PASS | report: sclet-ai-7b37cb3c recorded: sclet-ai-7b37cb3c |
| V10 | 报告中的 trade-off 判定与比较记录一致 | PASS | tradeoff: TRUE improved: TRUE |
| V11 | scVI 小节的数字 = scVI 运行记录里的数字 | PASS | mix: 0.989 pres: 0.06111 sec: 2877.99 |

总计：11/11 通过。
