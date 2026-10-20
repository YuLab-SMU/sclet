# M10 复验记录（CRC 恶性判定报告）

> 独立进程只读账本重算。时间 2026-10-20 21:28:50。
> 账本指纹 `sclet-ai-3c52eade`；对象 33694 × 1400。

| # | 断言 | 结果 | 细节 |
|---|---|---|---|
| V1 | 报告存在且含 4 条结论 | PASS | findings: 4 |
| V2 | 引用的证据 id 可解析 | PASS | refs: 3 |
| V3 | 每条记录等级 = 引擎重算等级 | PASS | 4/4 |
| V4 | CNV 节点记录了本次运行的 hmm/cutoff | PASS | hmm: TRUE cutoff: 1 |
| V5 | 四态计数 = colData 重算 | PASS | malignant: 0 non_malignant: 697 |
| V6 | 降级前 supported 保留且降级后肯定态为 0 | PASS | pre-demote supported: 4 post-demote affirmative: 0 |
| V7 | 签名阈值 = 全池中位数，且签名阳性数 = 证据节点 | PASS | cut=-0.00893 (median=TRUE) signature+=700/1400 |
| V8 | CNV+ 数 = 用记录阈值重算 | PASS | CNV+: 7 | threshold: 0.183 |
| V9 | 报告写出中位数判据、数据集级降级、CNV− ≠ 非恶性三条限制 | PASS | caveats present |
| V10 | RenderAIReport() 重渲染并通过自检 | PASS | self-check passed |

总计：10/10 通过。
