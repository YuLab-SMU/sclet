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

data_dir <- normalizePath(get_arg("data-dir", "regvelo_smoke_data"), mustWork = TRUE)
outdir <- normalizePath(get_arg("outdir", "regvelo_smoke_results"), mustWork = FALSE)
max_epochs <- as.integer(get_arg("max-epochs", "1"))
batch_size_arg <- get_arg("batch-size", NULL)
batch_size <- if (is.null(batch_size_arg) || !nzchar(batch_size_arg)) NULL else as.integer(batch_size_arg)
python <- get_arg("python", Sys.getenv("SCLET_REGVELO_PYTHON", unset = ""))
python <- if (nzchar(python)) python else NULL
sclet_path <- get_arg("sclet-path", Sys.getenv("SCLET_SOURCE", unset = ""))
sclet_path <- if (nzchar(sclet_path)) sclet_path else normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
seed <- as.integer(get_arg("seed", "1"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("RegVelo smoke test")
message("  data_dir: ", data_dir)
message("  outdir: ", outdir)
message("  sclet_path: ", sclet_path)
message("  max_epochs: ", max_epochs)

if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required.")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required to load sclet from source.")
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

torch_info <- tryCatch({
    torch <- reticulate::import("torch")
    list(
        version = torch$`__version__`,
        cuda_available = torch$cuda$is_available(),
        cuda_device_count = torch$cuda$device_count(),
        cuda_device_name = if (torch$cuda$is_available()) torch$cuda$get_device_name(0L) else NA_character_
    )
}, error = function(e) {
    list(error = conditionMessage(e))
})
message("Torch info:")
print(torch_info)
if (isTRUE(!is.null(torch_info$version))) {
    torch_major_minor <- sub("^([0-9]+\\.[0-9]+).*", "\\1", torch_info$version)
    if (utils::compareVersion(torch_major_minor, "2.6") >= 0) {
        warning(
            "RegVelo's published dependency metadata expects torch < 2.6. ",
            "This environment reports torch ", torch_info$version,
            ". The smoke test will continue, but a RegVelo import/runtime ",
            "failure here likely means this Python environment needs an older ",
            "PyTorch build."
        )
    }
}

module_versions <- list()
for (mod in c("anndata", "scanpy", "scvelo", "scvi", "regvelo")) {
    ok <- reticulate::py_module_available(mod)
    message("Python module ", mod, ": ", ok)
    if (!ok) {
        stop("Required Python module is missing: ", mod)
    }
    module_versions[[mod]] <- tryCatch({
        obj <- reticulate::import(mod)
        as.character(obj$`__version__`)
    }, error = function(e) NA_character_)
}
message("Python module versions:")
print(module_versions)

regvelo_r_path <- file.path(sclet_path, "R", "regvelo.R")
message("sclet R/regvelo.R: ", regvelo_r_path)
if (file.exists(regvelo_r_path)) {
    message("sclet R/regvelo.R mtime: ", file.info(regvelo_r_path)$mtime)
    regvelo_r_text <- readLines(regvelo_r_path, warn = FALSE)
    message(
        "sclet R/regvelo.R has debug marker: ",
        any(grepl("SCLET_REGVELO_DEBUG_VERSION=2026-07-01-velocity-keys", regvelo_r_text, fixed = TRUE))
    )
}

devtools::load_all(sclet_path, quiet = TRUE)

sce_path <- file.path(data_dir, "sce_regvelo_smoke.rds")
grn_path <- file.path(data_dir, "prior_grn.csv")
if (!file.exists(sce_path)) {
    stop("Missing SCE file: ", sce_path)
}
if (!file.exists(grn_path)) {
    stop("Missing GRN file: ", grn_path)
}

sce <- readRDS(sce_path)
grn <- utils::read.csv(grn_path, stringsAsFactors = FALSE)

message("Input dimensions: ", nrow(sce), " genes x ", ncol(sce), " cells")
message("GRN edges: ", nrow(grn))

set.seed(seed)
started <- Sys.time()
sce <- RunRegVelo(
    sce,
    grn = grn,
    reduction = "PCA",
    name = "regvelo_smoke",
    max_epochs = max_epochs,
    batch_size = batch_size,
    seed = seed,
    backend = "reticulate",
    python = python
)
finished <- Sys.time()

result_path <- file.path(outdir, "sce_regvelo_smoke_result.rds")
saveRDS(sce, result_path)

vel <- get_velocity(sce, id = "regvelo_smoke")
summary <- data.frame(
    key = c(
        "started",
        "finished",
        "elapsed_seconds",
        "n_cells",
        "n_genes",
        "method",
        "velocity_assay",
        "colData",
        "rowData",
        "torch_cuda_available",
        "torch_cuda_device_count",
        "torch_cuda_device_name",
        "python",
        paste0("py_", names(module_versions))
    ),
    value = c(
        as.character(started),
        as.character(finished),
        as.numeric(difftime(finished, started, units = "secs")),
        ncol(sce),
        nrow(sce),
        vel$method,
        vel$artifacts$velocity_assay,
        paste(vel$artifacts$colData, collapse = ";"),
        paste(vel$artifacts$rowData, collapse = ";"),
        torch_info$cuda_available %||% NA,
        torch_info$cuda_device_count %||% NA,
        torch_info$cuda_device_name %||% NA,
        py_config$python,
        unlist(module_versions, use.names = FALSE)
    )
)
utils::write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE)

message("Wrote result: ", result_path)
message("Wrote summary: ", file.path(outdir, "summary.csv"))
message("Done.")
