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
outdir <- normalizePath(get_arg("outdir", "regvelo_smoke_results/plots"), mustWork = FALSE)
sclet_path <- get_arg("sclet-path", Sys.getenv("SCLET_SOURCE", unset = ""))
sclet_path <- if (nzchar(sclet_path)) sclet_path else normalizePath(file.path(getwd()), mustWork = FALSE)
velocity_id <- get_arg("id", "regvelo_smoke")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required to load sclet from source.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to save plots.")
}

devtools::load_all(sclet_path, quiet = TRUE)

sce <- readRDS(result)
cell_magnitude <- VelocityMagnitude(sce, id = velocity_id, margin = "cell")
gene_magnitude <- VelocityMagnitude(sce, id = velocity_id, margin = "gene")

utils::write.csv(
    data.frame(cell = names(cell_magnitude), velocity_magnitude = cell_magnitude),
    file.path(outdir, "cell_velocity_magnitude.csv"),
    row.names = FALSE
)
utils::write.csv(
    data.frame(gene = names(gene_magnitude), velocity_magnitude = gene_magnitude),
    file.path(outdir, "gene_velocity_magnitude.csv"),
    row.names = FALSE
)

p <- plot_velocity_magnitude(sce, id = velocity_id)
ggplot2::ggsave(
    filename = file.path(outdir, "velocity_magnitude.png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 150
)

top_gene_plot <- plot_top_velocity_genes(sce, id = velocity_id, n = 20)
ggplot2::ggsave(
    filename = file.path(outdir, "top_velocity_genes.png"),
    plot = top_gene_plot,
    width = 6,
    height = 5,
    dpi = 150
)

message("Wrote: ", file.path(outdir, "cell_velocity_magnitude.csv"))
message("Wrote: ", file.path(outdir, "gene_velocity_magnitude.csv"))
message("Wrote: ", file.path(outdir, "velocity_magnitude.png"))
message("Wrote: ", file.path(outdir, "top_velocity_genes.png"))
