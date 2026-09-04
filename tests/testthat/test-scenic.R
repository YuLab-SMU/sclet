test_that("RunSCENIC validates inputs before reaching basilisk", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("basilisk")

    counts <- matrix(rpois(100, 5), nrow = 20)
    rownames(counts) <- paste0("gene", 1:20)
    colnames(counts) <- paste0("cell", 1:5)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    tfs <- tempfile(fileext = ".txt")
    writeLines(c("gene1", "gene2", "gene3"), tfs)
    motifs <- tempfile(fileext = ".tbl")
    writeLines(c("motif\tTF"), motifs)
    db <- tempfile(fileext = ".feather")
    writeLines("not a real feather", db)

    # Unknown assay
    expect_error(
        RunSCENIC(sce, tfs, motifs, db, assay_use = "counts2"),
        class = "sclet_missing_assay"
    )

    # Missing TFs file
    expect_error(
        RunSCENIC(sce, file.path(tempdir(), "nope_tfs.txt"), motifs, db),
        class = "sclet_file_not_found"
    )

    # Missing motif annotations file
    expect_error(
        RunSCENIC(sce, tfs, file.path(tempdir(), "nope_motifs.tbl"), db),
        class = "sclet_file_not_found"
    )

    # Missing database file
    expect_error(
        RunSCENIC(sce, tfs, motifs, file.path(tempdir(), "nope_db.feather")),
        class = "sclet_file_not_found"
    )

    # output_dir required for save_intermediate / resume
    expect_error(
        RunSCENIC(sce, tfs, motifs, db, save_intermediate = TRUE),
        class = "sclet_missing_output_dir"
    )
    expect_error(
        RunSCENIC(sce, tfs, motifs, db, resume = TRUE),
        class = "sclet_missing_output_dir"
    )
})

test_that("RunSCENIC stores the AUC altExp and registers a scenic record", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("basilisk")

    counts <- matrix(rpois(100, 5), nrow = 20)
    rownames(counts) <- paste0("gene", 1:20)
    colnames(counts) <- paste0("cell", 1:5)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    tfs <- tempfile(fileext = ".txt")
    writeLines(c("gene1", "gene2"), tfs)
    motifs <- tempfile(fileext = ".tbl")
    writeLines(c("motif\tTF"), motifs)
    db <- tempfile(fileext = ".feather")
    writeLines("not a real feather", db)

    outdir <- tempfile("scenic_")
    on.exit(unlink(outdir, recursive = TRUE), add = TRUE)

    captured <- NULL
    mock_run <- function(env, fun, ...) {
        captured <<- list(env = env, dots = list(...))
        # Return a cells x regulons data.frame as the Python callback would.
        data.frame(
            reg1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
            reg2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
            row.names = colnames(counts)
        )
    }

    out <- with_mocked_bindings(
        RunSCENIC(
            sce, tfs, motifs, db,
            assay_use = "counts",
            num_workers = 2L,
            name = "scenic_main",
            output_dir = outdir,
            save_intermediate = TRUE,
            resume = TRUE
        ),
        sclet_basilisk_run = mock_run,
        .package = "sclet"
    )

    # altExp stored with regulons (rows) x cells (cols)
    expect_true("SCENIC_AUC" %in% SingleCellExperiment::altExpNames(out))
    auc <- SingleCellExperiment::altExp(out, "SCENIC_AUC")
    expect_equal(nrow(auc), 2L)
    expect_equal(ncol(auc), 5L)

    # The basilisk callback received the checkpoint arguments.
    expect_equal(captured$dots$n_workers, 2L)
    expect_equal(captured$dots$out_dir, outdir)
    expect_true(captured$dots$save_mid)
    expect_true(captured$dots$resume)

    # The scenic state record is registered.
    rec <- get_scenic(out, id = "scenic_main")
    expect_equal(rec$params$num_workers, 2L)
    expect_true(rec$params$save_intermediate)
    expect_true(rec$params$resume)
    expect_equal(rec$artifacts$output_dir, outdir)
    expect_equal(rec$summary$n_regulons, 2L)
    expect_equal(rec$summary$n_cells, 5L)
})

test_that("sclet_scenic_env pins a full, tested SCENIC dependency set", {
    env <- sclet:::sclet_scenic_env
    pkgs <- env@packages
    expect_true("python=3.10" %in% pkgs)
    expect_true("pyscenic==0.12.1" %in% pkgs)
    expect_true("loompy==3.0.7" %in% pkgs)
    # Transitive pins needed for current pySCENIC compatibility (#27)
    expect_true("setuptools==80.9.0" %in% pkgs)
    expect_true("arboreto==0.1.6" %in% pkgs)
    expect_true("dask==2024.2.1" %in% pkgs)
    expect_true("distributed==2024.2.1" %in% pkgs)
    expect_true("dask-expr==0.5.3" %in% pkgs)
    expect_true("pyarrow==20.0.0" %in% pkgs)
})

test_that("sclet_env_setup_error raises an actionable sclet_env_setup_failed error", {
    env <- sclet:::sclet_scenic_env
    err <- tryCatch(
        sclet:::sclet_env_setup_error(
            env,
            simpleError("Error installing package(s): 'pip', 'wheel', 'setuptools'")
        ),
        error = function(e) e
    )
    expect_s3_class(err, "sclet_env_setup_failed")
    expect_match(conditionMessage(err), "sclet_scenic_env")
    expect_match(conditionMessage(err), "BASILISK_CUSTOM_PYTHON_3_10")
})

test_that("sclet_numpy_compat restores np.object when numpy is present", {
    skip_if_not(reticulate::py_available())
    has_numpy <- tryCatch({
        reticulate::import("numpy")
        TRUE
    }, error = function(e) FALSE)
    skip_if_not(has_numpy)

    # Restore the alias and ensure the shim is idempotent and error-free.
    sclet:::sclet_numpy_compat()
    np <- reticulate::import("numpy")
    expect_true(reticulate::py_has_attr(np, "object"))
})
