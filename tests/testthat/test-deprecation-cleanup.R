collect_warning_messages <- function(expr) {
    warnings <- character()
    withCallingHandlers(
        expr,
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    warnings
}

test_that("preprocess and batch helpers avoid legacy deprecation warnings", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("batchelor")

    set.seed(101)
    sce1 <- scuttle::mockSCE(ncells = 20, ngenes = 80)
    sce2 <- scuttle::mockSCE(ncells = 20, ngenes = 80)
    sce1$batch <- rep(c("a", "b"), each = 10)
    sce2$batch <- rep(c("c", "d"), each = 10)

    warning_messages <- collect_warning_messages({
        sce1 <- NormalizeData(sce1)
        sce1 <- FindVariableFeatures(sce1, nfeatures = 20, method = "scran")
        sce2 <- NormalizeData(sce2)
        sce2 <- FindVariableFeatures(sce2, nfeatures = 20, method = "scran")
        merged <- sce_merge(list(a = sce1, b = sce2))
        BatchRemover(merged, nHVG = 20, correct.all = FALSE)
    })

    deprecated_messages <- grep(
        "normalizeCounts|fitTrendVar|combineBlocks|getTopHVGs|deprecated|不再有用",
        warning_messages,
        value = TRUE,
        ignore.case = TRUE
    )

    expect_length(deprecated_messages, 0)
})
