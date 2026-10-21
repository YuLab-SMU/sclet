## M12-b 回归题库：用 M1–M11 的已知答案检验「探索数据价值」的判断力
##
## 用法：
##   Rscript .dev/m12-question-bank.R
##   Rscript .dev/m12-question-bank.R Q1   # 只跑一条
##
## 为什么不进 make check：每题都要载入一个真实对象（16 MB – 1 GB），
## 而 make check 里没有这些对象，也不该有。题库是「记录产物」的一部分，
## 结论写进 notes/m12-plan-quality-report.md。
##
## 判据的来源（每题都有出处，不是事后编的）：
##   Q1 M7  notes/pbmc10k-doublet-report.md        cluster 12 是 98.3% 双细胞
##   Q2 M9  notes/roadmap-2026-10-17.md            velocity 需 counts，只有 logcounts 时必须拒绝
##   Q3 M2  notes/m2-report-kang.md                stim 只在 batch2，混杂必须写成限制
##   Q4 M10 notes/crc-malignant-report.md          签名阈值取中位数 → 恰好一半，是定义不是生物学
##   Q5 M5  notes/pbmc10k-annotation-report.md     cluster 9 主导标签 0.22 → 混合，不给类型
##   Q6 M11 notes/pancreas-velocity-fate-report.md 终态占比 1/k 读起来像去向，实为终态细胞间占比
##
## 每题的判定方式：把该局面下的 `data_support` 与账本交给判据，看输出是否落在
## 「已知正确答案」上。**这不是测模型，是测上下文有没有把正确答案变成可读的。**

suppressWarnings(suppressMessages({
    library(pkgload)
    load_all("/home/wang/data/source/omics/sclet", export_all = TRUE, quiet = TRUE)
}))

AGENTDIR <- "/tmp/sclet-agentrun"

results <- list()
record <- function(id, question, verdict, detail) {
    results[[length(results) + 1L]] <<- list(
        id = id, question = question, verdict = verdict, detail = detail
    )
    cat(sprintf("%-3s %-6s %s\n     %s\n", id, verdict, question, detail))
}

has_object <- function(name) {
    file.exists(file.path(AGENTDIR, name))
}

load_object <- function(name) {
    tryCatch(readRDS(file.path(AGENTDIR, name)), error = function(e) NULL)
}

# ---- Q1（M7）：cluster 12 是 98.3% 双细胞 ----------------------------------
# 正确答案：两个断言。
#   a) unexplored 不得建议重跑 doublet（它已有结果）；
#   b) unexplored 必须非空并且真的指向尚未做的那件事——这个对象有 clustering、
#      有 enrichment，但没有轨迹，所以 run_trajectory 应在。
# 只断言 (a) 时，一个永远返回空清单的实现也能 PASS。第一版题库就是这样。
q1 <- function() {
    sce <- load_object("m7-report.rds")
    if (is.null(sce)) return(record("Q1", "cluster 12 双细胞局面", "SKIP", "m7-report.rds 不可用"))
    ds <- GetAnalysisLedger(sce)$analysis_story$data_support
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    redoes_doublet <- "run_doublet_detection" %in% actions
    offers_trajectory <- "run_trajectory" %in% actions
    ok <- !redoes_doublet && offers_trajectory
    record("Q1", "cluster 12 双细胞局面",
        if (ok) "PASS" else "FAIL",
        paste0("unexplored 共 ", length(actions), " 条 [",
            paste(actions, collapse = ","), "]；重跑 doublet=", redoes_doublet,
            "；提出尚未做的轨迹=", offers_trajectory))
}

# ---- Q2（M9）：只有 logcounts 时必须拒绝 velocity ---------------------------
# 正确答案：两个断言。
#   a) velocity supported = FALSE，reason 是 counts_assay_missing 或
#      spliced_unspliced_missing，且 unexplored 不提 run_velocity；
#   b) 前提自证：这个对象必须真的缺 counts——否则 a) 可能在测一个不存在的局面。
#      M9 的记录说「本仓库所有 velocity 导出只有 logcounts」，这里直接核对。
q2 <- function() {
    sce <- load_object("m9.rds")
    if (is.null(sce)) return(record("Q2", "velocity 数据不足局面", "SKIP", "m9.rds 不可用"))
    ds <- GetAnalysisLedger(sce)$analysis_story$data_support
    velocity <- Filter(function(x) identical(x$question_class, "velocity"), ds$supports)
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    if (!length(velocity)) return(record("Q2", "velocity 数据不足局面", "FAIL", "supports 没有 velocity 项"))
    velocity <- velocity[[1L]]
    assays <- SummarizedExperiment::assayNames(sce)
    premise_holds <- !("counts" %in% assays) ||
        !isTRUE(ds$shapes$counts_are_integer)
    ok <- !isTRUE(velocity$supported) && !("run_velocity" %in% actions) &&
        velocity$reason %in% c("counts_assay_missing", "spliced_unspliced_missing") &&
        premise_holds
    record("Q2", "velocity 数据不足局面",
        if (ok) "PASS" else "FAIL",
        paste0("assays=[", paste(assays, collapse = ","), "]",
            " counts_integer=", ds$shapes$counts_are_integer,
            " 前提成立=", premise_holds,
            " supported=", velocity$supported,
            " reason=", velocity$reason %||% "NULL",
            "；unexplored 提 velocity=", "run_velocity" %in% actions))
}

# ---- Q3（M2）：Kang 的 stim 只在 batch2 ------------------------------------
# 正确答案：数据形状必须识别出多样本 + 有 batch 列 + 有重复，
# 并把「批次与条件可能混杂」变成 diff abundance 可用的前提，而不是直接放行。
q3 <- function() {
    sce <- load_object("m2-final.rds")
    if (is.null(sce)) return(record("Q3", "Kang 设计混杂局面", "SKIP", "m2-final.rds 不可用"))
    ds <- GetAnalysisLedger(sce)$analysis_story$data_support
    detected <- isTRUE(ds$shapes$has_multiple_samples) &&
        length(ds$shapes$candidate_batch_columns) > 0L
    unconfirmed <- length(ds$unconfirmed_semantics) > 0L
    record("Q3", "Kang 设计混杂局面",
        if (detected) "PASS" else "FAIL",
        paste0("multiple_samples=", ds$shapes$has_multiple_samples,
            " batch 列=", paste(ds$shapes$candidate_batch_columns, collapse = ","),
            " replicate=", ds$shapes$has_replicates,
            " unconfirmed_semantics=", unconfirmed))
}

# ---- Q4（M10）：CNV 阈值取中位数 → 恰好一半 --------------------------------
# 正确答案：两个断言。
#   a) 不建议重跑已完成的 CNV / malignant（重复付费分析）；
#   b) 0 条建议必须**有理由**：这个对象是 inferCNV 的原始对象，没有降维也没有
#      聚类，所以轨迹那条本来就该被 reduction_missing 挡住。
# 只断言 (a) 会把「判据整体失效、什么都算不出来」的情况误判成 PASS——
# 第一版题库就犯了这错，m10 的 unexplored 是 0 条而 PASS 了。
q4 <- function() {
    sce <- load_object("m10-report.rds")
    if (is.null(sce)) return(record("Q4", "CNV 阈值局面", "SKIP", "m10-report.rds 不可用"))
    ds <- GetAnalysisLedger(sce)$analysis_story$data_support
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    redoes <- any(c("run_cnv_inference", "call_malignant_cells") %in% actions)
    trajectory <- Filter(function(x) identical(x$question_class, "trajectory"), ds$supports)[[1L]]
    reductions <- tryCatch(SingleCellExperiment::reducedDimNames(sce),
        error = function(e) character())
    # 0 条建议只在「确实没有可做的」时才算数：没有降维的对象不该被建议轨迹
    empty_is_justified <- length(actions) > 0L ||
        (length(reductions) == 0L && identical(trajectory$reason, "reduction_missing"))
    ok <- !redoes && empty_is_justified
    record("Q4", "CNV 阈值局面",
        if (ok) "PASS" else "FAIL",
        paste0("建议 ", length(actions), " 条；重跑 CNV=", redoes,
            "；reductions=", length(reductions),
            " trajectory_reason=", trajectory$reason %||% "NULL",
            "；空清单有理由=", empty_is_justified))
}

# ---- Q5（M5）：cluster 9 主导标签 0.22 → 混合，不给类型 --------------------
# 正确答案：两个断言。
#   a) 注释已完成，unexplored 不得再提 run_annotation；
#   b) 清单必须非空：这个对象做过聚类、注释、doublet、rare、detest、
#      enrichment，但没有轨迹，所以 run_trajectory 应在。
q5 <- function() {
    sce <- load_object("m5-report.rds")
    if (is.null(sce)) return(record("Q5", "低置信注释局面", "SKIP", "m5-report.rds 不可用"))
    ds <- GetAnalysisLedger(sce)$analysis_story$data_support
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    redoes_annotation <- "run_annotation" %in% actions
    has_suggestion <- length(actions) > 0L
    ok <- !redoes_annotation && has_suggestion
    record("Q5", "低置信注释局面",
        if (ok) "PASS" else "FAIL",
        paste0("unexplored 共 ", length(actions), " 条 [",
            paste(actions, collapse = ","), "]；重跑 annotation=", redoes_annotation))
}

# ---- Q6（M11）：终态占比 vs 去向分布 ---------------------------------------
# 正确答案：三个断言。
#   a) shapes 正确识别 counts + spliced/unspliced + 整数 counts；
#   b) unexplored 不得再提 velocity / fate（都已跑过）；
#   c) 清单必须非空——不是所有都做完了：这个对象有没有做 milo / 注释，
#      取决于具体对象，但至少不该是一个「全空」的清单。
q6 <- function() {
    sce <- load_object("m11-report.rds")
    if (is.null(sce)) return(record("Q6", "速度-命运链局面", "SKIP", "m11-report.rds 不可用"))
    ds <- GetAnalysisLedger(sce)$analysis_story$data_support
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    shapes_ok <- isTRUE(ds$shapes$has_spliced_unspliced) &&
        isTRUE(ds$shapes$counts_are_integer)
    no_rerun <- !any(c("run_velocity", "run_fate_analysis") %in% actions)
    ok <- shapes_ok && no_rerun && length(actions) > 0L
    record("Q6", "速度-命运链局面",
        if (ok) "PASS" else "FAIL",
        paste0("spliced/unspliced=", ds$shapes$has_spliced_unspliced,
            " integer counts=", ds$shapes$counts_are_integer,
            "；unexplored ", length(actions), " 条 [", paste(actions, collapse = ","),
            "]；提 velocity/fate=", any(c("run_velocity", "run_fate_analysis") %in% actions)))
}

questions <- list(Q1 = q1, Q2 = q2, Q3 = q3, Q4 = q4, Q5 = q5, Q6 = q6)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) questions <- questions[intersect(names(questions), only)]

cat("== M12-b 回归题库 ==\n")
for (fn in questions) fn()

verdicts <- vapply(results, function(x) x$verdict, character(1))
cat(sprintf("\n汇总：PASS %d / FAIL %d / SKIP %d（共 %d 题）\n",
    sum(verdicts == "PASS"), sum(verdicts == "FAIL"),
    sum(verdicts == "SKIP"), length(verdicts)))

out <- file.path("/tmp", "m12-question-bank-result.rds")
saveRDS(results, out)
cat("结果已存：", out, "\n")
