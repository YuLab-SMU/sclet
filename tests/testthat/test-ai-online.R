test_that("real aisdk integration is opt-in and key-gated", {
    skip_if_not_installed("aisdk")
    skip_if(
        !nzchar(Sys.getenv("DEEPSEEK_API_KEY")),
        "DEEPSEEK_API_KEY is not configured"
    )
    skip_if(
        !identical(Sys.getenv("SCLET_RUN_ONLINE_TESTS"), "true"),
        "set SCLET_RUN_ONLINE_TESTS=true to enable online provider tests"
    )

    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(c(1, 0, 3, 2, 0, 1, 4, 1, 0, 2, 1, 3), nrow = 4, ncol = 3))
    )
    result <- AIStatus(sce, model = "deepseek:deepseek-chat")
    expect_s3_class(result, "sclet_ai_result")
    expect_true(isTRUE(result$metadata$structured_output))
})
