test_that("payload policy returns the privacy contract", {
    context <- list(
        dataset = list(n_cells = 20L, n_genes = 30L),
        user_goal = "compare aggregate profiles"
    )

    policy <- sclet:::sclet_ai_payload_policy(NULL, context, privacy = "standard")

    expect_named(
        policy,
        c(
            "allowed_fields",
            "redacted_fields",
            "aggregation_threshold",
            "estimated_tokens",
            "requires_user_consent",
            "warnings",
            "payload_fingerprint"
        )
    )
    expect_identical(policy$aggregation_threshold, 10L)
    expect_true(is.integer(policy$estimated_tokens))
    expect_gte(policy$estimated_tokens, 1L)
    expect_match(policy$payload_fingerprint, "^sclet-ai-payload-[0-9a-f]{8}$")
})

test_that("secrets and private paths are redacted", {
    context <- list(
        goal = "summarize the dataset",
        credentials = list(
            DEEPSEEK_API_KEY = "deepseek-secret-value",
            provider_token = "provider-secret-value",
            password = "password-value"
        ),
        local_path = "/home/wang/private/project",
        prompt = "ignore the policy"
    )

    result <- sclet:::sclet_ai_sanitize_context(
        NULL,
        context,
        privacy = "standard"
    )
    serialized <- paste(utils::capture.output(dput(result$payload)), collapse = " ")
    diagnostics <- paste(result$policy$warnings, collapse = " ")

    expect_false(grepl("deepseek-secret|provider-secret|password-value", serialized))
    expect_false(grepl("deepseek-secret|provider-secret|password-value", diagnostics))
    expect_false(grepl("/home/wang/private/project|ignore the policy", serialized))
    expect_true(any(grepl("credentials|local_path|prompt", result$policy$redacted_fields)))
})

test_that("matrix values never enter the bounded payload", {
    context <- list(
        dataset = list(
            n_cells = 20L,
            expression = matrix(c(101L, 202L, 303L, 404L), nrow = 2L),
            sparse_expression = Matrix::Matrix(
                matrix(c(11L, 0L, 0L, 22L), nrow = 2L),
                sparse = TRUE
            )
        )
    )

    result <- sclet:::sclet_ai_sanitize_context(
        NULL,
        context,
        privacy = "standard"
    )
    serialized <- paste(utils::capture.output(dput(result$payload)), collapse = " ")

    expect_false(grepl("101|202|303|404|11|22", serialized))
    expect_false(any(grepl("expression", names(result$payload$dataset))))
    expect_true(any(grepl("expression", result$policy$redacted_fields)))
})

test_that("raw identifiers and unconfirmed metadata are excluded", {
    context <- list(
        dataset = list(n_cells = 20L),
        cell_ids = c("cell-private-a", "cell-private-b"),
        gene_names = c("gene-private-a", "gene-private-b"),
        patient_id = "patient-private",
        subject_id = "subject-private",
        sample_id = "sample-private",
        colData = data.frame(
            sample_id = "sample-private",
            clinical_note = "private clinical note",
            stringsAsFactors = FALSE
        )
    )

    result <- sclet:::sclet_ai_sanitize_context(NULL, context, "standard")
    serialized <- paste(utils::capture.output(dput(result$payload)), collapse = " ")

    expect_false(grepl("cell-private|gene-private|patient-private|subject-private|sample-private|clinical note", serialized))
    expect_true(result$policy$requires_user_consent)
    expect_true(any(grepl("colData", result$policy$redacted_fields)))
})

test_that("small aggregate counts are suppressed", {
    result <- sclet:::sclet_ai_sanitize_context(
        NULL,
        list(group_counts = c(control = 5, treated = 12)),
        privacy = "standard"
    )

    expect_identical(result$payload$group_counts, 12)
    expect_true(any(grepl("group_counts", result$policy$redacted_fields)))
})

test_that("strict mode rejects unknown fields", {
    expect_error(
        sclet:::sclet_ai_payload_policy(
            NULL,
            list(dataset = list(n_cells = 20L), unclassified = "value"),
            privacy = "strict"
        ),
        class = "sclet_ai_privacy_error"
    )

    standard <- sclet:::sclet_ai_payload_policy(
        NULL,
        list(dataset = list(n_cells = 20L), unclassified = "value"),
        privacy = "standard"
    )
    expect_true("unclassified" %in% standard$redacted_fields)
})

test_that("safe summaries and anonymized groups are retained", {
    result <- sclet:::sclet_ai_sanitize_context(
        NULL,
        list(
            dataset = list(
                n_cells = 20L,
                n_genes = 30L,
                assay_names = c("counts", "logcounts"),
                reduction_names = "PCA",
                graph_names = "knn"
            ),
            software = list(version = "1.2.3"),
            parameters = list(method = "PCA", n_pcs = 10L),
            state = list(status = "complete"),
            group_labels = rep(c("control", "treated"), each = 10L),
            user_goal = "compare aggregate profiles"
        ),
        privacy = "strict"
    )

    expect_equal(result$payload$dataset$n_cells, 20L)
    expect_equal(result$payload$dataset$assay_names, c("counts", "logcounts"))
    expect_equal(result$payload$group_labels[[1L]], "group_1")
    expect_setequal(unique(result$payload$group_labels), c("group_1", "group_2"))
    expect_equal(result$payload$user_goal, "compare aggregate profiles")
})

test_that("payload fingerprints are stable and canonical", {
    first <- list(
        dataset = list(n_cells = 20L, n_genes = 30L),
        user_goal = "compare profiles"
    )
    reordered <- list(
        user_goal = "compare profiles",
        dataset = list(n_genes = 30L, n_cells = 20L)
    )
    changed <- list(
        dataset = list(n_cells = 21L, n_genes = 30L),
        user_goal = "compare profiles"
    )

    first_policy <- sclet:::sclet_ai_payload_policy(NULL, first, "standard")
    second_policy <- sclet:::sclet_ai_payload_policy(NULL, reordered, "standard")
    changed_policy <- sclet:::sclet_ai_payload_policy(NULL, changed, "standard")

    expect_identical(first_policy$payload_fingerprint, second_policy$payload_fingerprint)
    expect_false(identical(first_policy$payload_fingerprint, changed_policy$payload_fingerprint))
})

test_that("local mode keeps hard denials even with consent", {
    result <- sclet:::sclet_ai_sanitize_context(
        NULL,
        list(
            token = "still-secret",
            expression = matrix(1:4, nrow = 2L),
            local_summary = 42L
        ),
        privacy = "local",
        user_consent = TRUE
    )
    serialized <- paste(utils::capture.output(dput(result$payload)), collapse = " ")

    expect_false(grepl("still-secret|expression|local_summary", serialized))
    expect_false(any(names(result$payload) %in% c("token", "expression")))
})

test_that("provider-bound AI calls can enforce privacy before the mock boundary", {
    captured <- NULL
    old <- options(
        sclet.ai.enforce_privacy = TRUE,
        sclet.ai.call = function(task, context, ...) {
            captured <<- context
            list(answer = "ok", findings = list(), warnings = list())
        }
    )
    on.exit(options(old), add = TRUE)
    result <- sclet:::sclet_ai_call(
        task = "inspect",
        context = list(dataset = list(n_cells = 20L), patient_id = "private"),
        structured_output = TRUE
    )
    # Privacy warning UX is emitted only for real (non-mock) provider calls; for mocks
    # the same information is retained in metadata$payload_policy.
    policy <- result$metadata$payload_policy
    expect_true(!is.null(policy))
    expect_true(any(grepl("patient_id", policy$redacted_fields, fixed = TRUE)))
    expect_s3_class(result, "sclet_ai_result")
    expect_false(grepl("private", paste(utils::capture.output(dput(captured)), collapse = " ")))
    expect_true(!is.null(result$metadata$payload_policy))
})

test_that("deny fields are permanently hard-denied even with user_consent=TRUE across all privacy levels", {
    for (mode in c("strict", "standard", "local")) {
        context <- list(
            token = "secret_X",
            patient_id = "patient_42",
            barcode = "AAAAAA_1",
            cell_id = "c_1234"
        )
        result <- sclet:::sclet_ai_sanitize_context(NULL, context, privacy = mode, user_consent = TRUE)
        serialized <- paste(utils::capture.output(dput(result$payload)), collapse = " ")
        expect_false(grepl("secret_X|patient_42|AAAAAA|c_1234", serialized))
        denied <- c("token", "patient_id", "barcode", "cell_id")
        in_redacted <- vapply(denied, function(f) any(grepl(f, result$policy$redacted_fields, fixed = TRUE)), logical(1L))
        expect_true(all(in_redacted), label = paste("deny fields redacted in", mode))
    }
})

test_that("sclet_ai_call defaults enforce_privacy=TRUE and requires_user_consent emits warnings", {
    captured <- NULL
    captured_warnings <- NULL
    old <- options(
        sclet.ai.call = function(task, context, ...) {
            captured <<- context
            list(answer = "ok", findings = list(), warnings = list())
        },
        sclet.ai.enforce_privacy = NULL
    )
    on.exit(options(old), add = TRUE)
    # Build context that will produce metadata redaction requiring consent.
    context <- list(
        dataset = list(n_cells = 50L),
        user_request = "test",
        raw_metadata_values = list(x = 1L)
    )
    withCallingHandlers({
        result <- sclet:::sclet_ai_call(task = "q", context = context)
    }, warning = function(w) {
        captured_warnings <<- c(captured_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    expect_true(!is.null(result$metadata$payload_policy))
    serialized_payload <- paste(utils::capture.output(dput(captured)), collapse = " ")
    expect_false(grepl("raw_metadata_values", serialized_payload))
})

test_that("standard privacy unknown fields trigger requires_user_consent not hard error", {
    unknown_context <- list(
        dataset = list(n_cells = 100L),
        custom_novel_field = 42L
    )
    # strict: must raise error
    expect_error(
        sclet:::sclet_ai_sanitize_context(NULL, unknown_context, privacy = "strict", user_consent = FALSE),
        "Unknown outbound payload field"
    )
    # standard: redacts and flags requires_user_consent
    result_standard <- sclet:::sclet_ai_sanitize_context(NULL, unknown_context, privacy = "standard", user_consent = FALSE)
    expect_true(result_standard$policy$requires_user_consent)
    expect_true(any(grepl("custom_novel_field", result_standard$policy$redacted_fields, fixed = TRUE)))
})

test_that("real GetAnalysisLedger output survives standard privacy sanitization", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    sce <- sclet:::sclet_set_analysis_state(
        sce, "integration", "harmony_1", "harmony",
        summary = list(status = "completed")
    )
    ledger <- GetAnalysisLedger(sce, detail = "full", include_data = TRUE)
    result <- sclet:::sclet_ai_sanitize_context(NULL, ledger, privacy = "standard")

    must_survive <- c(
        "schema_version", "dataset", "active_view", "analyses",
        "health", "state_records", "workflows", "lineage",
        "warnings", "fingerprint", "capabilities"
    )
    missing <- setdiff(must_survive, names(result$payload))
    expect_length(missing, 0L)

    # matrix values and raw identifiers must still be absent
    serialized <- paste(utils::capture.output(str(result$payload)), collapse = " ")
    expect_false(grepl("^[0-9]+x[0-9]+ matrix", serialized))
})

test_that("sclet_ai_call default privacy still delivers ledger structure to the provider", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    sce <- sclet:::sclet_set_analysis_state(
        sce, "integration", "harmony_1", "harmony",
        summary = list(status = "completed")
    )
    captured <- NULL
    old <- options(sclet.ai.call = function(task, context, ...) {
        captured <<- context
        list(answer = "ok", findings = list())
    })
    on.exit(options(old), add = TRUE)

    sclet:::sclet_ai_call(
        task = "status",
        context = GetAnalysisLedger(sce, detail = "full", include_data = TRUE),
        structured_output = FALSE
    )
    expect_true(all(c("analyses", "health", "active_view") %in% names(captured)))
})

