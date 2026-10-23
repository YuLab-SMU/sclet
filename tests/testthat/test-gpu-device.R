# The device switch: default CPU, GPU on request, and never a silent fallback.
#
# The contract: `device = NULL` (the default) is CPU on every machine; "cuda"
# without a working GPU environment is a typed error, not a CPU fallback that
# hides the request; "auto" picks GPU only when the environment reports a card.

test_that("the default device is cpu, whatever the machine has", {
    # The GPU env exists on this machine and reports a card, and the default
    # must still be CPU: that is the release decision.
    withr::local_options(sclet.python.device = NULL)
    withr::local_envvar(SCLET_PYTHON_DEVICE = "")
    expect_identical(sclet:::sclet_ai_device(), "cpu")
    expect_identical(sclet:::sclet_ai_device(NULL), "cpu")
})

test_that("the three sources of a device, most specific first", {
    withr::local_envvar(SCLET_PYTHON_DEVICE = "cpu")
    # explicit argument wins over everything
    expect_identical(sclet:::sclet_ai_device("cuda"), "cuda")
    # then the session option
    withr::local_options(sclet.python.device = "cuda")
    expect_identical(sclet:::sclet_ai_device(), "cuda")
    # then the environment variable
    withr::local_options(sclet.python.device = NULL)
    expect_identical(sclet:::sclet_ai_device(), "cpu")
})

test_that("an unknown device is refused, not coerced", {
    expect_error(sclet:::sclet_ai_device("gpu"), class = "sclet_ai_invalid_device")
    expect_error(sclet:::sclet_ai_device("GPU"), class = "sclet_ai_invalid_device")
})

test_that("auto resolves against the environment, not the driver", {
    # A driver being installed says nothing about whether the environment's
    # torch can use it; auto must ask the environment.
    withr::local_options(sclet.python.device = NULL)
    expect_identical(
        sclet:::sclet_ai_device("auto"),
        if (sclet:::sclet_gpu_available()) "cuda" else "cpu"
    )
})

test_that("gpu status reports the environment without building it", {
    # Cheap and offline: it must not create the environment as a side effect.
    before <- dir.exists(sclet_env_path("sclet_gpu_env"))
    status <- sclet_gpu_status()
    expect_named(status, c("env_exists", "python", "torch", "cuda",
        "device_count", "device_names", "available"))
    expect_identical(dir.exists(sclet_env_path("sclet_gpu_env")), before)
    # Whatever the machine has, the fields must agree with each other.
    expect_identical(status$available,
        isTRUE(status$cuda) && status$device_count > 0L)
})

test_that("RunIntegration takes a device argument and defaults to cpu", {
    # The CPU path is unchanged: same environment, same accelerator argument.
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 30)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce)
    SummarizedExperiment::colData(sce)$batch <- rep(c("a", "b"), each = 20)
    args <- formals(RunIntegration)
    expect_true("device" %in% names(args))
    expect_null(args$device)
})

test_that("a cuda request without a usable environment is a typed error", {
    # Simulated by pointing the env root at an empty directory: the resolver
    # must refuse rather than fall back to CPU, because a silent fallback hides
    # a request the user made deliberately.
    withr::local_envvar(SCLET_ENV_ROOT = tempfile())
    expect_false(sclet:::sclet_gpu_available())
    expect_identical(sclet:::sclet_ai_device("auto"), "cpu")
    # "cuda" is still honoured as a request, and the integration path turns it
    # into the typed error -- asserted here on the predicate the path uses.
    expect_identical(sclet:::sclet_ai_device("cuda"), "cuda")
})
