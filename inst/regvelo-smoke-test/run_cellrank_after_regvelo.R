#!/usr/bin/env Rscript

sclet_smoke_original_libpaths <- .libPaths()
sclet_smoke_rlib <- if (nzchar(Sys.getenv("HOME"))) {
    file.path(Sys.getenv("HOME"), "Rlib")
} else {
    path.expand("~/Rlib")
}
dir.create(sclet_smoke_rlib, showWarnings = FALSE, recursive = TRUE)
dir.create(path.expand("~/Rlib"), showWarnings = FALSE, recursive = TRUE)
.libPaths('~/Rlib')
.libPaths(unique(c(sclet_smoke_rlib, path.expand("~/Rlib"), sclet_smoke_original_libpaths, .libPaths())))

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
get_arg <- function(name, default = NULL) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (!length(hit)) {
        return(default)
    }
    sub(paste0("^--", name, "="), "", hit[[1]])
}

result <- normalizePath(
    get_arg("result", "regvelo_smoke_results/sce_regvelo_smoke_result.rds"),
    mustWork = TRUE
)
outdir <- normalizePath(get_arg("outdir", "regvelo_smoke_cellrank"), mustWork = FALSE)
sclet_path <- get_arg("sclet-path", Sys.getenv("SCLET_SOURCE", unset = ""))
sclet_path <- if (nzchar(sclet_path)) sclet_path else normalizePath(file.path(getwd()), mustWork = FALSE)
sclet_path <- normalizePath(sclet_path, mustWork = TRUE)
python <- get_arg("python", Sys.getenv("SCLET_REGVELO_PYTHON", unset = ""))
python <- if (nzchar(python)) python else NULL
velocity_id <- get_arg("velocity-id", "regvelo_smoke")
trajectory_id <- get_arg("trajectory-id", "cellrank_regvelo_smoke")
reduction <- get_arg("reduction", "PCA")
cluster_key <- get_arg("cluster-key", "")
cluster_key <- if (nzchar(cluster_key)) cluster_key else NULL
k <- as.integer(get_arg("k", "4"))
seed <- as.integer(get_arg("seed", "1"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("RegVelo -> CellRank smoke test")
message("  result: ", result)
message("  outdir: ", outdir)
message("  sclet_path: ", sclet_path)
message("  velocity_id: ", velocity_id)
message("  trajectory_id: ", trajectory_id)
message("  reduction: ", reduction)

if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required.")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required to load sclet from source.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to save plots.")
}

cellrank_source <- file.path(sclet_path, "R", "cellrank.R")
message("sclet R/cellrank.R: ", cellrank_source)
if (!file.exists(cellrank_source)) {
    stop("Cannot find sclet source file: ", cellrank_source)
}
message("sclet R/cellrank.R mtime: ", file.info(cellrank_source)$mtime)
cellrank_source_text <- readLines(cellrank_source, warn = FALSE)
required_source_markers <- c(
    "sclet_prepare_cellrank_velocity_input",
    "backend = c(\"basilisk\", \"reticulate\")",
    "velocity_assay = velocity_assay"
)
missing_source_markers <- required_source_markers[
    !vapply(required_source_markers, function(marker) {
        any(grepl(marker, cellrank_source_text, fixed = TRUE))
    }, logical(1))
]
if (length(missing_source_markers)) {
    stop(
        "The sclet source at ", sclet_path,
        " does not contain the RegVelo-aware CellRank changes. Missing marker(s): ",
        paste(missing_source_markers, collapse = ", "),
        ". Upload or pull the updated source tree before rerunning this job."
    )
}

if (!is.null(python)) {
    message("Binding reticulate to Python/conda env: ", python)
    if (file.exists(python)) {
        reticulate::use_python(python, required = TRUE)
    } else {
        reticulate::use_condaenv(python, required = TRUE)
    }
}

py_config <- reticulate::py_config()
message("Python: ", py_config$python)
message("Python version: ", py_config$version)
for (mod in c("anndata", "scanpy", "scvelo", "cellrank")) {
    ok <- reticulate::py_module_available(mod)
    message("Python module ", mod, ": ", ok)
    if (!ok) {
        stop("Required Python module is missing: ", mod)
    }
}

devtools::load_all(sclet_path, quiet = TRUE)
run_cellrank_formals <- names(formals(get("RunCellRank", envir = asNamespace("sclet"))))
missing_run_cellrank_args <- setdiff(c("velocity_id", "backend", "python"), run_cellrank_formals)
has_cellrank_helper <- exists(
    "sclet_prepare_cellrank_velocity_input",
    envir = asNamespace("sclet"),
    inherits = FALSE
)
if (length(missing_run_cellrank_args) || !has_cellrank_helper) {
    stop(
        "Loaded sclet::RunCellRank() is not the RegVelo-aware version. ",
        "Loaded formals: ", paste(run_cellrank_formals, collapse = ", "), ". ",
        "sclet_path: ", sclet_path, ". ",
        "R/cellrank.R mtime: ", file.info(cellrank_source)$mtime, ". ",
        "Make sure the server copy of sclet includes the latest R/cellrank.R."
    )
}

sce <- readRDS(result)
if (!velocity_id %in% names(sclet:::sclet_get_state_records(sce, "velocity"))) {
    stop("Velocity id not found in object: ", velocity_id)
}
velocity_record <- get_velocity(sce, id = velocity_id)
velocity_assay <- velocity_record$artifacts$velocity_assay
message("Velocity assay recorded for ", velocity_id, ": ", velocity_assay %||% "<none>")
if (is.null(velocity_assay) || !velocity_assay %in% SummarizedExperiment::assayNames(sce)) {
    stop(
        "Velocity id ", velocity_id, " does not point to a stored velocity assay in the SCE. ",
        "Available assays: ", paste(SummarizedExperiment::assayNames(sce), collapse = ", ")
    )
}
if (!reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
    stop("Reduction not found in object: ", reduction)
}

if (is.null(cluster_key)) {
    cluster_key <- "regvelo_smoke_cluster"
    emb <- as.matrix(SingleCellExperiment::reducedDim(sce, reduction))
    set.seed(seed)
    k <- max(2L, min(k, nrow(emb)))
    SummarizedExperiment::colData(sce)[[cluster_key]] <- paste0("cluster", stats::kmeans(emb, centers = k)$cluster)
    message("Created temporary cluster column: ", cluster_key, " (k=", k, ")")
} else if (!cluster_key %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("cluster_key not found in colData: ", cluster_key)
}

started <- Sys.time()
sce <- RunCellRank(
    sce,
    reduction = reduction,
    cluster_key = cluster_key,
    velocity_id = velocity_id,
    backend = "reticulate",
    python = python,
    name = trajectory_id
)
finished <- Sys.time()

result_path <- file.path(outdir, "sce_regvelo_cellrank_result.rds")
saveRDS(sce, result_path)

cellrank <- get_cellrank(sce, id = trajectory_id)
summary <- data.frame(
    key = c(
        "started",
        "finished",
        "elapsed_seconds",
        "n_cells",
        "n_genes",
        "velocity_id",
        "trajectory_id",
        "cluster_key",
        "reduction",
        "backend",
        "python",
        "terminal_state_col",
        "fate_probability_cols",
        "has_lineage_drivers"
    ),
    value = c(
        as.character(started),
        as.character(finished),
        as.numeric(difftime(finished, started, units = "secs")),
        ncol(sce),
        nrow(sce),
        velocity_id,
        trajectory_id,
        cluster_key,
        reduction,
        cellrank$inputs$backend %||% "reticulate",
        py_config$python,
        cellrank$artifacts$terminal_state_col %||% "",
        paste(cellrank$artifacts$fate_probability_cols %||% character(), collapse = ";"),
        cellrank$artifacts$has_lineage_drivers %||% FALSE
    )
)
utils::write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE)

cellrank_diagnostics <- CellRankSummary(sce, id = trajectory_id)
utils::write.csv(
    cellrank_diagnostics$terminal_states,
    file.path(outdir, "cellrank_terminal_state_counts.csv"),
    row.names = FALSE
)
utils::write.csv(
    cellrank_diagnostics$fate_probabilities,
    file.path(outdir, "cellrank_fate_probability_summary.csv"),
    row.names = FALSE
)
velocity_fate_correlation <- VelocityFateCorrelation(
    sce,
    velocity_id = velocity_id,
    trajectory_id = trajectory_id
)
utils::write.csv(
    velocity_fate_correlation,
    file.path(outdir, "velocity_fate_correlation.csv"),
    row.names = FALSE
)

terminal_plot <- plot_fate_terminal_states(sce, id = trajectory_id, reduction = reduction)
ggplot2::ggsave(
    filename = file.path(outdir, "cellrank_terminal_states.png"),
    plot = terminal_plot,
    width = 6,
    height = 5,
    dpi = 150
)

fate_cols <- cellrank$artifacts$fate_probability_cols %||% character()
if (length(fate_cols)) {
    fate_plot <- plot_fate_probability(sce, id = trajectory_id, fate = 1, reduction = reduction)
    ggplot2::ggsave(
        filename = file.path(outdir, "cellrank_fate_probability_1.png"),
        plot = fate_plot,
        width = 6,
        height = 5,
        dpi = 150
    )
    velocity_fate_plot <- plot_velocity_fate_correlation(
        sce,
        fate = 1,
        velocity_id = velocity_id,
        trajectory_id = trajectory_id
    )
    ggplot2::ggsave(
        filename = file.path(outdir, "velocity_fate_correlation.png"),
        plot = velocity_fate_plot,
        width = 6,
        height = 5,
        dpi = 150
    )
}

message("Wrote result: ", result_path)
message("Wrote summary: ", file.path(outdir, "summary.csv"))
message("Wrote diagnostics: ", file.path(outdir, "cellrank_terminal_state_counts.csv"))
message("Wrote diagnostics: ", file.path(outdir, "cellrank_fate_probability_summary.csv"))
message("Wrote diagnostics: ", file.path(outdir, "velocity_fate_correlation.csv"))
message("Wrote plot: ", file.path(outdir, "cellrank_terminal_states.png"))
if (length(fate_cols)) {
    message("Wrote plot: ", file.path(outdir, "cellrank_fate_probability_1.png"))
    message("Wrote plot: ", file.path(outdir, "velocity_fate_correlation.png"))
}
message("Done.")
