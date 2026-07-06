# Prepare a small RegVelo smoke-test dataset.
#
# Preferred usage from an interactive R terminal:
#
#   source("inst/regvelo-smoke-test/prepare_data.R")
#   prepare_regvelo_smoke_data(outdir = "regvelo_smoke_data")
#
# For backwards compatibility this file can still be called by Rscript with
# --outdir=..., --n-cells=..., --n-genes=..., --seed=..., --source-h5ad=...

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

prepare_regvelo_smoke_data <- function(
    outdir = "regvelo_smoke_data",
    n_cells = 300L,
    n_genes = 800L,
    seed = 1L,
    source_h5ad = NULL
) {
    outdir <- normalizePath(outdir, mustWork = FALSE)
    n_cells <- as.integer(n_cells)
    n_genes <- as.integer(n_genes)
    seed <- as.integer(seed)

    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    set.seed(seed)

    message("Preparing RegVelo smoke-test data in: ", outdir)

    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'SingleCellExperiment' is required.")
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Package 'SummarizedExperiment' is required.")
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required.")
    }

    if (!is.null(source_h5ad) && nzchar(source_h5ad)) {
        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop("Package 'reticulate' is required when source_h5ad is used.")
        }
        if (grepl("^/path/to/", source_h5ad) || !file.exists(source_h5ad)) {
            stop(
                "source_h5ad does not point to an existing file: ", source_h5ad, "\n",
                "Create or download an h5ad first. For example:\n",
                "  python inst/regvelo-smoke-test/download_smoke_h5ad.py ",
                "--out regvelo_smoke_source/pancreas.h5ad\n",
                "Then call prepare_regvelo_smoke_data(source_h5ad = ",
                "'regvelo_smoke_source/pancreas.h5ad', ...)."
            )
        }
        message("Reading source h5ad through Python anndata: ", source_h5ad)
        ad <- reticulate::import("anndata", convert = FALSE)
        adata <- ad$read_h5ad(source_h5ad)
        adata <- sclet_smoke_subset_adata(adata, n_cells = n_cells, n_genes = n_genes)
        sce <- sclet_smoke_adata_to_sce(adata)
    } else {
        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop("Package 'reticulate' is required to download scvelo sample data.")
        }
        message("Downloading scvelo pancreas sample data via Python...")
        scv <- reticulate::import("scvelo", convert = FALSE)
        adata <- scv$datasets$pancreas()
        adata <- sclet_smoke_subset_adata(adata, n_cells = n_cells * 2L, n_genes = n_genes * 2L)
        sce <- sclet_smoke_adata_to_sce(adata)
    }

    assay_names <- SummarizedExperiment::assayNames(sce)
    if (!all(c("spliced", "unspliced") %in% assay_names)) {
        stop(
            "The input data must contain assays named 'spliced' and 'unspliced'. ",
            "Available assays: ", paste(assay_names, collapse = ", ")
        )
    }

    keep_genes <- rownames(sce)
    if (is.null(keep_genes)) {
        keep_genes <- paste0("gene", seq_len(nrow(sce)))
        rownames(sce) <- keep_genes
    }

    spliced <- SummarizedExperiment::assay(sce, "spliced")
    gene_score <- Matrix::rowSums(Matrix::Matrix(spliced, sparse = TRUE))
    gene_score[is.na(gene_score)] <- 0
    gene_order <- order(gene_score, decreasing = TRUE)
    gene_order <- gene_order[seq_len(min(n_genes, length(gene_order)))]

    cell_score <- Matrix::colSums(Matrix::Matrix(spliced, sparse = TRUE))
    cell_score[is.na(cell_score)] <- 0
    cell_order <- order(cell_score, decreasing = TRUE)
    cell_order <- cell_order[seq_len(min(n_cells, length(cell_order)))]

    sce <- sce[gene_order, cell_order]

    if (!"PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
        if (!requireNamespace("scater", quietly = TRUE)) {
            stop("Package 'scater' is required to compute PCA for the smoke-test subset.")
        }
        if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
            if (!requireNamespace("scuttle", quietly = TRUE)) {
                stop("Package 'scuttle' is required to create logcounts.")
            }
            sce <- scuttle::logNormCounts(sce, assay.type = "spliced")
        }
        sce <- scater::runPCA(sce, exprs_values = "logcounts", ncomponents = 20)
    }

    genes <- rownames(sce)
    n_regulators <- min(80L, max(10L, floor(length(genes) * 0.1)))
    regulators <- genes[seq_len(n_regulators)]
    targets <- genes

    edge_count <- min(2000L, length(targets) * length(regulators))
    grn <- data.frame(
        target = sample(targets, edge_count, replace = TRUE),
        regulator = sample(regulators, edge_count, replace = TRUE),
        weight = runif(edge_count, min = 0.05, max = 1),
        stringsAsFactors = FALSE
    )
    grn <- grn[grn$target != grn$regulator, , drop = FALSE]
    grn <- grn[!duplicated(paste(grn$target, grn$regulator, sep = "\r")), , drop = FALSE]

    sce_path <- file.path(outdir, "sce_regvelo_smoke.rds")
    grn_path <- file.path(outdir, "prior_grn.csv")
    manifest_path <- file.path(outdir, "manifest.csv")

    saveRDS(sce, sce_path)
    utils::write.csv(grn, grn_path, row.names = FALSE)

    manifest <- data.frame(
        key = c("n_cells", "n_genes", "n_regulators", "n_edges", "sce", "grn"),
        value = c(
            ncol(sce),
            nrow(sce),
            length(unique(grn$regulator)),
            nrow(grn),
            sce_path,
            grn_path
        )
    )
    utils::write.csv(manifest, manifest_path, row.names = FALSE)

    message("Wrote: ", sce_path)
    message("Wrote: ", grn_path)
    message("Done.")

    invisible(list(
        sce = sce,
        grn = grn,
        manifest = manifest,
        files = list(sce = sce_path, grn = grn_path, manifest = manifest_path)
    ))
}

sclet_smoke_subset_adata <- function(adata, n_cells, n_genes) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required.")
    }

    reticulate::py_run_string(
        paste(
            "def sclet_smoke_subset_adata(adata, n_cells, n_genes):",
            "    n_cells = min(int(n_cells), adata.n_obs)",
            "    n_genes = min(int(n_genes), adata.n_vars)",
            "    return adata[:n_cells, :n_genes].copy()",
            sep = "\n"
        )
    )
    reticulate::py$sclet_smoke_subset_adata(
        adata,
        as.integer(n_cells),
        as.integer(n_genes)
    )
}

sclet_smoke_adata_to_sce <- function(adata) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required.")
    }
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'SingleCellExperiment' is required.")
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Package 'SummarizedExperiment' is required.")
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required.")
    }

    py_attr <- function(x, name, default = NULL) {
        if (reticulate::py_has_attr(x, name)) {
            reticulate::py_get_attr(x, name)
        } else {
            default
        }
    }
    py_get_item <- function(x, key, default = NULL) {
        tryCatch(reticulate::py_get_item(x, key), error = function(e) default)
    }
    py_to_r <- function(x) reticulate::py_to_r(x)

    reticulate::py_run_string(
        paste(
            "def sclet_smoke_axis_names(adata):",
            "    return [adata.obs_names.astype(str).tolist(), adata.var_names.astype(str).tolist()]",
            sep = "\n"
        )
    )
    axis_names <- reticulate::py$sclet_smoke_axis_names(adata)
    obs_names <- as.character(py_to_r(axis_names[[1]]))
    var_names <- as.character(py_to_r(axis_names[[2]]))
    layers <- py_attr(adata, "layers")

    as_gene_cell_matrix <- function(x, layer_name) {
        if (is.null(x)) {
            stop("AnnData layer '", layer_name, "' was not found.")
        }
        if (reticulate::py_has_attr(x, "toarray")) {
            x <- reticulate::py_call(py_attr(x, "toarray"))
        }
        mat <- as.matrix(py_to_r(x))
        if (length(obs_names) == 1L && length(var_names) == 1L && length(dim(mat)) == 2L) {
            obs_names <<- paste0("cell", seq_len(nrow(mat)))
            var_names <<- paste0("gene", seq_len(ncol(mat)))
        }
        expected_cell_gene <- c(length(obs_names), length(var_names))
        expected_gene_cell <- rev(expected_cell_gene)
        if (identical(dim(mat), expected_cell_gene)) {
            mat <- t(mat)
        } else if (!identical(dim(mat), expected_gene_cell)) {
            stop(
                "Unexpected AnnData layer '", layer_name, "' dimensions: got ",
                paste(dim(mat), collapse = " x "),
                ", expected cells x genes (", paste(expected_cell_gene, collapse = " x "),
                ") or genes x cells (", paste(expected_gene_cell, collapse = " x "), ")."
            )
        }
        rownames(mat) <- var_names
        colnames(mat) <- obs_names
        Matrix::Matrix(mat, sparse = TRUE)
    }

    spliced <- as_gene_cell_matrix(py_get_item(layers, "spliced"), "spliced")
    unspliced <- as_gene_cell_matrix(py_get_item(layers, "unspliced"), "unspliced")
    SingleCellExperiment::SingleCellExperiment(
        assays = list(spliced = spliced, unspliced = unspliced)
    )
}

sclet_prepare_data_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
    get_arg <- function(name, default = NULL) {
        hit <- grep(paste0("^--", name, "="), args, value = TRUE)
        if (!length(hit)) {
            return(default)
        }
        sub(paste0("^--", name, "="), "", hit[[1]])
    }

    list(
        outdir = get_arg("outdir", "regvelo_smoke_data"),
        n_cells = as.integer(get_arg("n-cells", "300")),
        n_genes = as.integer(get_arg("n-genes", "800")),
        seed = as.integer(get_arg("seed", "1")),
        source_h5ad = get_arg("source-h5ad", NULL)
    )
}

sclet_prepare_data_called_as_script <- function() {
    any(grepl("(^|[/\\\\])prepare_data[.]R$", sub("^--file=", "", commandArgs(FALSE))))
}

if (sclet_prepare_data_called_as_script() && identical(Sys.getenv("SCLET_PREPARE_DATA_NO_RUN"), "")) {
    do.call(prepare_regvelo_smoke_data, sclet_prepare_data_cli_args())
}
