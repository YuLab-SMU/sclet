# ai_report_crc_malignant

- Recorded: 2026-10-20 21:28:22.034798
- Data: 33694 genes x 1400 cells
- Ledger fingerprint at record time: `sclet-ai-066927c6`; now: `sclet-ai-3c52eade`
- Findings: 4

## Summary

CRC 队列的恶性判定完成：CNV 阈值 0.183，CNV+ 7 个；四态最终 malignant 0、likely_malignant 0、undetermined 703、non_malignant 697。

## Findings (claim levels derived by the engine)

| Finding | Evidence cited | Recorded level | Engine audit |
|---|---|---|---|
| inferCNV 在 1400 个细胞上得到 burden（中位 0.0971，P90 0.129），参考细胞 300 个，阈值取自参考 burden（P99+3×MAD = 0.183），越线者为 CNV+：7 个。 | `ev:cnv_cnv`, `ev:malignant_malignant` | associated | ok |
| 四态调用：malignant 0、likely_malignant 0、undetermined 703、non_malignant 697；签名分阈值 -0.0089250542844088（中位数默认：1），签名阳性 700 个细胞。 | `ev:malignant_malignant_identity`, `ev:cnv_cnv` | associated | ok |
| 降级前的判定保留在 malignant_call_pre_demote：supported 4、probable 3；降级后这两类都变成 undetermined。 | `ev:malignant_malignant_identity` | observed | ok |
| 这些细胞里哪些是恶性，需要在这个队列上验证 CNV 线的灵敏度；本次没有独立的恶金标准。 | - | hypothesis | ok |

## Evidence

| Evidence | claim_level | parents | dependency_group |
|---|---|---|---|
| `ev:cnv_cnv` | associated | - | dc0bf290ae75 |
| `ev:malignant_malignant` | associated | - | dc0bf290ae75 |
| `ev:malignant_malignant_identity` | associated | - | dc0bf290ae75 |

## Claim ceiling

- Ceiling: `associated` (1 independent line(s) of support)
- Reason: all supporting evidence resolves to 1 independent line(s) of support; a claim cannot exceed 'associated'

## Warnings and recommendations

- Warning: identity 竞争检验在本次数据集上触发了降级：降级前 supported 4 个，降级后 0 个。
- Warning: 签名阈值默认取全池中位数，因此「签名阳性」是相对陈述。
- Warning: CNV 阈值由参考细胞推导（P99+3×MAD），在本队列上只放出少量 CNV+。
- Recommendation: 若要比较不同阈值，固定 hmm 与 cutoff 后重跑，并把两个比例一起报告。
- Recommendation: 若要判断 CNV 线的灵敏度，需要 DNA 级或病理级仲裁。

## Limitations stated in the findings

- **inferCNV 在 1400 个细胞上得到 burden（中位 0.0971，P90 0.129），参考细胞 300 个，阈值取自参考 burden（P99+3×MAD = 0.183），越线者为 CNV+：7 个。**: 本次运行 parameter_hmm = TRUE，cutoff = 1；HMM 关掉时 burden 来自非 HMM profile，是另一种测量而不是同一种的加速版。; CNV 阴性不等于非恶性：阈值保守时「看不到 CNV」会被读成否定证据。
- **四态调用：malignant 0、likely_malignant 0、undetermined 703、non_malignant 697；签名分阈值 -0.0089250542844088（中位数默认：1），签名阳性 700 个细胞。**: 签名阈值取全池中位数时，签名阳性按定义接近一半细胞，比例本身不是生物学信号。; identity 竞争检验是数据集级的：一个分离统计量会把所有 supported/probable 降级。
- **降级前的判定保留在 malignant_call_pre_demote：supported 4、probable 3；降级后这两类都变成 undetermined。**: 这是对象状态的直接读数，不是生物学结论。
- **这些细胞里哪些是恶性，需要在这个队列上验证 CNV 线的灵敏度；本次没有独立的恶金标准。**: 没有 DNA 级仲裁数据；队列的 CNV 稀少程度未知。

## Recorded analyses

### annotation `cnv`

- Method: inferCNV per-cell copy-number burden

### annotation `malignant`

- Method: four-state malignant-cell call from CNV + signature evidence

### annotation `malignant_identity`

- Method: four-state malignant-cell call from CNV + signature evidence

## Executed steps (from the execution records)

### ai_execution_plan_m10_cnv_20261020201516_121z44jk_20261020211704

| step | action | status | seconds |
|---|---|---|---|
| c1 | run_cnv_inference | completed | 3707.1 |

### ai_execution_plan_m10_malignant_20261020212728_lak72jng_20261020212729

| step | action | status | seconds |
|---|---|---|---|
| m1 | call_malignant_cells | completed | 0.2 |

### ai_execution_plan_m10_malignant_identity_20261020212735_u0qcmg3v_20261020212737

| step | action | status | seconds |
|---|---|---|---|
| m2 | call_malignant_cells | completed | 0.55 |


## Verification

- Citations resolve: PASS of 3
- Claim levels re-derived: every recorded level matches the engine
- Result: all checks passed

> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model.
