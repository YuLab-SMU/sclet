# Phenotype evidence: a registry of independent lines that can support (or
# refuse) a per-cell phenotype call.
#
# The malignant-cell work in R/malignancy-cnv.R and R/malignancy-combine.R made
# the point concrete: a per-cell phenotype label is worth exactly as much as the
# number of INDEPENDENT observations behind it, and the most common failure is
# one observation dressed up as a verdict. Malignancy happens to have an
# unusually good second line (copy number is orthogonal to any transcriptional
# identity programme). Most phenotypes do not, and a method that pretends
# otherwise is worse than one that says so.
#
# So the registry does not decide anything. Each provider reports, per cell,
# `positive` / `negative` / `unavailable` for ONE line of evidence, records what
# it measured, and states how independent it is from the signature that defined
# the phenotype. CallPhenotypeCells() then counts lines: two agreeing lines
# support a positive call, one line makes it `probable`, no line or a conflict
# leaves it `undetermined`.
#
# A provider is a plain list so a user can register their own without subclassing
# anything:
#
#   list(
#     name     = "cnv",
#     label    = "Copy-number burden",
#     decide   = function(object, ...) list(state = c("positive", "negative", ...)),
#     available = function(object) list(available = TRUE, issues = character()),
#     dependence = "independent"   # or "shared_with_signature"
#   )
#
# `dependence` is the field that keeps the claim honest: a line derived from the
# same expression matrix as the signature is NOT a second observation, and the
# caller must say so. Counting it as independent would inflate the call, which is
# exactly the error the iterative method in R/malignancy.R makes when it re-derives
# its own signature.

#' Registry of phenotype evidence providers
#'
#' Holds the evidence providers [CallPhenotypeCells()] can draw on, keyed by name.
#' Built-in providers are registered at load time; a caller can add their own.
#'
#' @param action One of `"list"`, `"get"`, `"register"`, `"unregister"`.
#' @param name Provider name.
#' @param provider A provider list, see \sQuote{Providers}.
#' @param include Optional character vector restricting the returned names.
#'
#'
#' @return `action = "list"` returns the provider names; `"get"` returns one
#'   provider; `"register"` and `"unregister"` return `NULL` invisibly.
#'
#' @section Providers:
#' A provider is a list with:
#' \describe{
#'   \item{`name`}{Unique key.}
#'   \item{`label`}{Human-readable description used in messages.}
#'   \item{`available`}{`function(object, ...)` returning
#'     `list(available = TRUE/FALSE, issues = character())`. A provider that
#'     cannot run must say `FALSE` with a reason; it must never report
#'     `available = TRUE` and then decide `negative` for "nothing measured".}
#'   \item{`decide`}{`function(object, ...)` returning
#'     `list(state = character_vector_of_cells, detail = list(...))` where each
#'     state is `"positive"`, `"negative"` or `"unavailable"`.}
#'   \item{`dependence`}{`"independent"` or `"shared_with_signature"`. Only
#'     independent lines count toward a two-line call.}
#'   \item{`source`}{Optional string naming the MEASUREMENT a provider reads.
#'     Providers declaring the same `source` are one line of evidence, not
#'     several: two wrappers over one measurement are not corroborating
#'     observations. Use it whenever you register more than one provider over the
#'     same data.}
#' }
#'
#' @seealso [CallPhenotypeCells()], [RunCNVInference()]
#' @export
#' @examples
#' PhenotypeEvidenceProviders("list")
#' PhenotypeEvidenceProviders("get", name = "cnv")$label
PhenotypeEvidenceProviders <- function(action = c("list", "get", "register", "unregister"),
                                       name = NULL,
                                       provider = NULL,
                                       include = NULL) {
    action <- match.arg(action)
    registry <- sclet_phenotype_registry()
    switch(
        action,
        list = {
            nms <- names(registry)
            if (!is.null(include)) {
                unknown <- setdiff(include, nms)
                if (length(unknown) > 0) {
                    stop("Unknown evidence provider(s): ", paste(unknown, collapse = ", "),
                         ". Available: ", paste(nms, collapse = ", "), ".")
                }
                nms <- intersect(nms, include)
            }
            nms
        },
        get = {
            if (is.null(name) || !name %in% names(registry)) {
                stop("No evidence provider named '", name %||% "", "'.")
            }
            registry[[name]]
        },
        register = {
            if (is.null(provider) || !is.list(provider)) {
                stop("`provider` must be a provider list.")
            }
            check_phenotype_provider(provider)
            reg <- sclet_phenotype_registry()
            # A silent overwrite would let two different definitions share one
            # name, so replacing an existing provider must be asked for.
            if (provider$name %in% names(reg)) {
                stop(
                    "An evidence provider named '", provider$name, "' is already ",
                    "registered. Unregister it first if you mean to replace it."
                )
            }
            reg[[provider$name]] <- provider
            sclet_phenotype_registry_set(reg)
            invisible(NULL)
        },
        unregister = {
            if (is.null(name)) {
                stop("`name` is required for unregister.")
            }
            reg <- sclet_phenotype_registry()
            reg[[name]] <- NULL
            sclet_phenotype_registry_set(reg)
            invisible(NULL)
        }
    )
}

#' Check a phenotype evidence provider
#'
#' Validates a provider list before it is registered. Returns `TRUE` invisibly on
#' success and names every problem it finds otherwise.
#'
#' @param provider A provider list, see [PhenotypeEvidenceProviders()].
#'
#' @return `TRUE` invisibly, or an error naming the defects.
#' @export
#' @examples
#' provider <- list(
#'     name = "my_line",
#'     label = "An independent measurement",
#'     dependence = "independent",
#'     available = function(object, ...) list(available = TRUE, issues = character()),
#'     decide = function(object, ...) list(state = rep("positive", ncol(object)))
#' )
#' check_phenotype_provider(provider)
check_phenotype_provider <- function(provider) {
    problems <- character(0)
    if (!is.list(provider)) {
        stop("`provider` must be a list.")
    }
    for (field in c("name", "label", "available", "decide")) {
        if (is.null(provider[[field]])) {
            problems <- c(problems, paste0("missing field `", field, "`"))
        }
    }
    if (!is.null(provider$name) && (!is.character(provider$name) || !nzchar(provider$name))) {
        problems <- c(problems, "`name` must be a non-empty string")
    }
    for (fn in c("available", "decide")) {
        if (!is.null(provider[[fn]]) && !is.function(provider[[fn]])) {
            problems <- c(problems, paste0("`", fn, "` must be a function"))
        }
    }
    dep <- provider$dependence %||% "independent"
    if (!dep %in% c("independent", "shared_with_signature")) {
        problems <- c(problems, paste0(
            "`dependence` must be \"independent\" or \"shared_with_signature\", not \"",
            dep, "\""
        ))
    }
    if (!is.null(provider$source) &&
        (!is.character(provider$source) || length(provider$source) != 1 || !nzchar(provider$source))) {
        problems <- c(problems, "`source` must be a single non-empty string when supplied")
    }
    if (length(problems) > 0) {
        stop("Invalid phenotype evidence provider: ", paste(problems, collapse = "; "), ".")
    }
    invisible(TRUE)
}

#' Call a phenotype from independent lines of evidence
#'
#' Reports, for every cell, how well a phenotype call is supported and by which
#' lines -- rather than returning a single label that cannot say how sure it is.
#'
#' @param object A `SingleCellExperiment` object.
#' @param signature Character vector of genes, a named list with a `tumor`/
#'   `malignant` element, or the output of `FindMalignantSeedGenes()`. The
#'   signature defines the phenotype; it is scored and used as one line.
#' @param signature_score Optional `colData` column already holding the score.
#' @param evidence Character vector of provider names to draw on. Defaults to
#'   every registered provider; providers that are unavailable are reported, not
#'   silently dropped.
#' @param transcriptomic_cut Optional numeric cutoff on the signature score.
#'   `NULL` uses the median. The choice is recorded in the result either way.
#' @param competing Optional named list of competing signatures (`list(
#'    hypothesis_name = genes, ...)`). For each, the call reports how well the
#'    signature separates the evidence-positive cells from the rest, which is
#'    the question "is this the axis we think it is?". The groups are defined
#'    by the EVIDENCE (at least one line positive vs none), never by the
#'    resulting call: comparing against the call would make the test circular,
#'    since the signature selects the cells it is then asked to explain.
#' @param competing_demote Numeric. Separation on the competing signature
#'    above which affirmative calls (`supported`, `probable`) are demoted to
#'    `undetermined`. A `refuted` call is never demoted: it was never
#'    evidence-positive, so the axis question does not apply to it.
#' @param assay Assay for scoring. Defaults to `"logcounts"`.
#' @param name Analysis record id. Defaults to `"phenotype_call"`.
#' @param ... Arguments passed to each provider's `decide()`.
#'
#'
#' @return `object` with `colData` columns `<name>_call`, `<name>_claim`,
#'   `<name>_n_lines`, `<name>_signature_score` and one `<name>_<provider>_state`
#'   per provider, plus a state record of type `annotation` under `name`. The
#'   full result, including per-provider availability and issues, is at
#'   `metadata(object)$phenotype_call[[name]]`.
#'
#' @section The call and the claim:
#' `_call` is one of `supported`, `probable`, `undetermined`, `refuted`:
#' \describe{
#'   \item{`supported`}{Two or more lines are positive (the signature counts as
#'     one).}
#'   \item{`probable`}{Exactly one line is positive.}
#'   \item{`undetermined`}{A line says yes and another says no, no line is
#'     usable, or a competing signature explains the split.}
#'   \item{`refuted`}{No line says yes and at least one says no.}
#' }
#'
#' The asymmetry between a missing line and a negative one is deliberate. A low
#' signature score means *this line did not fire*, not *the phenotype is absent*;
#' only a provider that positively reports `negative` counts as a no. That is why
#' a signature-negative cell with one positive provider is `probable` rather than
#' `undetermined`, and why `refuted` requires an explicit negative.
#'
#' `_claim` reuses the package's claim ladder (`observed`, `measured`,
#' `associated`, `consistent_with`) rather than inventing a second vocabulary for
#' confidence: a supported call is `consistent_with` the phenotype, a probable
#' call is `associated` with it, and nothing here reaches `observed` because a
#' per-cell label is never a direct observation of cell state.
#'
#' @section What this does not do:
#' It does not iterate. The iterative search in [ClassifyMalignantEpithelium()] is
#' appropriate when the cell population itself is unknown; a call built on named
#' evidence is not improved by re-deriving its own signature, and doing so is the
#' positive feedback the reference method is criticised for. Phenotypes whose
#' evidence providers are all `shared_with_signature` cannot reach `supported` at
#' all -- that is deliberate, not an oversight.
#'
#' @seealso [PhenotypeEvidenceProviders()], [CallMalignantCells()]
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' sce <- NormalizeData(sce)
#' sce <- CallPhenotypeCells(sce, signature = rownames(sce)[1:10])
#' table(SummarizedExperiment::colData(sce)$phenotype_call_call)
CallPhenotypeCells <- function(object,
                               signature = NULL,
                               signature_score = NULL,
                               evidence = NULL,
                               transcriptomic_cut = NULL,
                               competing = NULL,
                               competing_demote = 0.1,
                               assay = NULL,
                               name = "phenotype_call",
                               ...) {
    if (!nzchar(name)) {
        stop("`name` must be a non-empty prefix.")
    }
    cdata <- SummarizedExperiment::colData(object)
    prefix <- paste0(name, "_")

    # ---- signature line ------------------------------------------------------
    if (!is.null(signature_score)) {
        if (!signature_score %in% colnames(cdata)) {
            stop("`signature_score` column '", signature_score, "' not found.")
        }
        sig_score <- as.numeric(cdata[[signature_score]])
        signature_used <- list(status = "supplied", column = signature_score)
    } else {
        genes <- sclet_phenotype_resolve_signature(signature)
        if (is.null(assay)) {
            assay <- if ("logcounts" %in% SummarizedExperiment::assayNames(object)) "logcounts" else NULL
        }
        if (is.null(assay)) {
            stop("No `logcounts` assay; pass `assay` or provide `signature_score`.")
        }
        expr <- sclet_malignancy_dense(SummarizedExperiment::assay(object, assay))
        bins <- sclet_gene_expression_bins(expr, n_bins = 25)
        sig_score <- sclet_add_module_score(expr, genes, bins = bins, n_ctrl = 100,
                                            exclude = genes)
        signature_used <- list(status = "scored", assay = assay, genes = length(genes))
    }
    cut_used <- if (is.null(transcriptomic_cut)) {
        list(value = stats::median(sig_score, na.rm = TRUE), source = "median")
    } else {
        list(value = transcriptomic_cut, source = "caller")
    }
    sig_positive <- !is.na(sig_score) & sig_score >= cut_used$value

    # ---- evidence lines ------------------------------------------------------
    available_providers <- PhenotypeEvidenceProviders("list")
    wanted <- if (is.null(evidence)) available_providers else evidence
    unknown <- setdiff(wanted, available_providers)
    if (length(unknown) > 0) {
        stop("Unknown evidence provider(s): ", paste(unknown, collapse = ", "),
             ". Available: ", paste(available_providers, collapse = ", "), ".")
    }

    provider_states <- list()
    provider_report <- list()
    for (pname in wanted) {
        provider <- PhenotypeEvidenceProviders("get", name = pname)
        avail <- tryCatch(
            provider$available(object, ...),
            error = function(e) list(available = FALSE, issues = conditionMessage(e))
        )
        if (!isTRUE(avail$available)) {
            provider_states[[pname]] <- rep("unavailable", nrow(cdata))
            provider_report[[pname]] <- list(
                state = "unavailable",
                label = provider$label,
                dependence = provider$dependence %||% "independent",
                source = provider$source,
                issues = avail$issues %||% character(0)
            )
            next
        }
        decided <- tryCatch(
            provider$decide(object, ...),
            error = function(e) list(state = NULL, error = conditionMessage(e))
        )
        if (is.null(decided$state)) {
            provider_states[[pname]] <- rep("unavailable", nrow(cdata))
            provider_report[[pname]] <- list(
                state = "unavailable",
                label = provider$label,
                dependence = provider$dependence %||% "independent",
                source = provider$source,
                issues = decided$error %||% "provider returned no state"
            )
            next
        }
        st <- as.character(decided$state)
        if (length(st) != nrow(cdata)) {
            stop("Provider '", pname, "' returned ", length(st),
                 " states for ", nrow(cdata), " cells.")
        }
        # A provider that names its states by cell is telling us the order; use
        # that instead of assuming positional. Without this, a provider that
        # returns its states in a different order than colData silently assigns
        # every cell the wrong verdict.
        if (!is.null(names(decided$state))) {
            nm <- names(decided$state)
            if (!setequal(nm, rownames(cdata))) {
                stop("Provider '", pname, "' returned named states whose names do not ",
                     "match colData cells (", length(setdiff(nm, rownames(cdata))),
                     " unknown, ", length(setdiff(rownames(cdata), nm)), " missing).")
            }
            st <- unname(st[match(rownames(cdata), nm)])
        }
        bad <- setdiff(st, c("positive", "negative", "unavailable"))
        if (length(bad) > 0) {
            stop("Provider '", pname, "' returned invalid state(s): ",
                 paste(unique(bad), collapse = ", "), ".")
        }
        provider_states[[pname]] <- st
        provider_report[[pname]] <- list(
            state = "measured",
            label = provider$label,
            dependence = provider$dependence %||% "independent",
            source = provider$source,
            issues = character(0),
            detail = decided$detail %||% NULL
        )
    }

    # ---- line accounting -----------------------------------------------------
    #
    # Providers that declare the same `source` are ONE line of evidence, not
    # several: two wrappers over one measurement are not corroborating
    # observations. Providers with no `source` are their own line.
    #
    # This closes the hole where registering the same thing twice under two names
    # (or two providers that silently read the same matrix) manufactured a
    # two-line call. The registry cannot see inside a provider, so `source` is a
    # declaration the provider makes; what it CAN enforce is that a declaration
    # is honoured, and that `dependence` still gates whether a line counts.
    independent <- names(Filter(
        function(p) identical(p$dependence %||% "independent", "independent") &&
            identical(p$state, "measured"),
        provider_report
    ))
    line_ids <- vapply(independent, function(p) {
        src <- provider_report[[p]]$source
        if (is.null(src) || !nzchar(src)) p else paste0("src:", src)
    }, character(1))
    unique_lines <- unique(line_ids)

    count_line <- function(target) {
        if (length(unique_lines) == 0) {
            return(rep(0L, nrow(cdata)))
        }
        Reduce(`+`, lapply(unique_lines, function(ln) {
            members <- independent[line_ids == ln]
            # A line is positive if ANY of its members is positive...
            pos <- Reduce(`|`, lapply(members, function(p) provider_states[[p]] == "positive"),
                          init = rep(FALSE, nrow(cdata)))
            # ...and negative only if it is positive nowhere and some member says
            # no. Agreement inside a line is not dissent.
            neg <- !pos & Reduce(`|`, lapply(members, function(p) provider_states[[p]] == "negative"),
                                 init = rep(FALSE, nrow(cdata)))
            as.integer(pos == target)
        }), init = rep(0L, nrow(cdata)))
    }
    n_positive <- count_line(TRUE)
    n_negative <- count_line(FALSE)

    # The signature is itself one line of evidence, and the only one that is
    # always available. Count it alongside the independent providers: a positive
    # signature plus a positive CNV is two agreeing lines, which is the whole
    # point of the two-line call.
    sig_line <- as.integer(sig_positive)
    n_yes <- sig_line + n_positive
    n_no <- n_negative

    call <- rep("undetermined", nrow(cdata))
    call[n_yes >= 2] <- "supported"
    call[n_yes == 1] <- "probable"
    call[n_yes == 0 & n_no > 0] <- "refuted"
    # A conflict -- at least one line says yes and at least one says no -- must
    # not be broken by counting more yes than no.
    conflicted <- n_yes > 0 & n_no > 0
    call[conflicted] <- "undetermined"
    # No usable provider and a negative signature is still refuted by that one
    # observation; with nothing at all it stays undetermined.
    call[conflicted & n_no == 0] <- "undetermined"

    total_lines <- n_yes + n_no

    # ---- competing hypotheses ------------------------------------------------
    #
    # The comparison groups come from the EVIDENCE-POSITIVE cells, not from
    # `call`. Using `call` would make this circular: the signature selects the
    # cells, then the competing programme is asked whether it explains "the cells
    # the signature selected", which it must, because the signature chose them.
    # The axis question is therefore asked of the evidence alone:
    #
    #   in_scope : at least one line is positive
    #   residual : every line is negative, or nothing was measured
    #
    # A competing programme separating those two groups says the evidence is
    # picking up that programme; it does not say so about cells that were already
    # refused, which must keep their refusal.
    call_pre_demote <- call          # defined even when no competing signature runs
    competing_report <- NULL
    if (!is.null(competing)) {
        if (is.null(assay)) {
            assay <- if ("logcounts" %in% SummarizedExperiment::assayNames(object)) "logcounts" else NULL
        }
        if (!is.null(assay)) {
            expr <- sclet_malignancy_dense(SummarizedExperiment::assay(object, assay))
            bins <- sclet_gene_expression_bins(expr, n_bins = 25)
            in_scope <- n_yes > 0
            residual <- !in_scope
            competing_report <- list()
            for (hname in names(competing)) {
                genes <- intersect(competing[[hname]], rownames(object))
                if (length(genes) == 0) {
                    competing_report[[hname]] <- list(status = "no_genes_in_object")
                    next
                }
                cscore <- sclet_add_module_score(expr, genes, bins = bins, n_ctrl = 100,
                                                 exclude = genes)
                usable <- !is.na(cscore)
                sep <- if (sum(in_scope & usable) > 2 && sum(residual & usable) > 2) {
                    mean(cscore[in_scope & usable]) - mean(cscore[residual & usable])
                } else {
                    NA_real_
                }
                competing_report[[hname]] <- list(
                    status = "measured",
                    genes = length(genes),
                    separation = sep,
                    n_in_scope = sum(in_scope & usable),
                    n_residual = sum(residual & usable),
                    demote = !is.na(sep) && abs(sep) > competing_demote
                )
            }
            demoting <- names(Filter(function(x) isTRUE(x$demote), competing_report))
            call_pre_demote <- call
            if (length(demoting) > 0) {
                # Only the affirmative states are demoted. `refuted` is a refusal
                # the competing programme cannot overturn -- it was never
                # evidence-positive, so the axis question does not apply to it.
                #
                # The demote is DATASET-LEVEL: one separation statistic removes
                # every affirmative cell, so on a dataset where the competing
                # programme separates the two groups the malignant state becomes
                # unreachable (measured by a user on GSE150290, 2026-10-20: 29
                # supported cells all demoted, final malignant count 0). The
                # pre-demote call is therefore kept rather than overwritten --
                # a reader must be able to see what the evidence said before the
                # axis check, and the report has to show both.
                call[call %in% c("supported", "probable")] <- "undetermined"
            }
        }
    }

    claim <- ifelse(call == "supported", "consistent_with",
                    ifelse(call == "probable", "associated", "hypothesis"))

    # ---- write back ----------------------------------------------------------
    cdata[[paste0(prefix, "call_pre_demote")]] <- factor(call_pre_demote,
        levels = c("supported", "probable", "undetermined", "refuted"))
    cdata[[paste0(prefix, "call")]] <- factor(call, levels = c("supported", "probable",
                                                              "undetermined", "refuted"))
    cdata[[paste0(prefix, "claim")]] <- claim
    cdata[[paste0(prefix, "n_lines")]] <- as.integer(total_lines)
    cdata[[paste0(prefix, "signature_score")]] <- sig_score
    for (pname in names(provider_states)) {
        cdata[[paste0(prefix, pname, "_state")]] <- provider_states[[pname]]
    }
    SummarizedExperiment::colData(object) <- cdata

    result <- list(
        call = table(factor(call, levels = c("supported", "probable", "undetermined", "refuted"))),
        n_cells = nrow(cdata),
        signature = signature_used,
        signature_cut = cut_used,
        n_signature_positive = as.integer(sum(sig_positive)),
        n_pre_demote_supported = as.integer(sum(call_pre_demote == "supported")),
        n_pre_demote_probable = as.integer(sum(call_pre_demote == "probable")),
        providers = provider_report,
        independent_lines = independent,
        competing = competing_report,
        thresholds = list(transcriptomic_cut = cut_used$value,
                          competing_demote = competing_demote),
        name = name
    )
    # Keyed by run id, mirroring metadata$phenotype_call[[name]], so two calls do
    # not overwrite each other and a reader can address one by name.
    existing <- sclet_get_analysis(object, "phenotype_call")
    if (is.null(existing)) existing <- list()
    existing[[name]] <- result
    object <- sclet_set_analysis(object, "phenotype_call", existing)
    md <- S4Vectors::metadata(object)
    md$phenotype_call <- if (is.null(md$phenotype_call)) list() else md$phenotype_call
    md$phenotype_call[[name]] <- result
    S4Vectors::metadata(object) <- md

    object <- sclet_set_state_record(
        object = object,
        type = "annotation",
        id = name,
        value = list(
            method = "phenotype call from independent lines of evidence",
            inputs = list(
                signature = signature_used$status,
                providers = names(provider_report)
            ),
            artifacts = list(columns = paste0(prefix, c("call", "claim", "n_lines"))),
            params = list(transcriptomic_cut = cut_used$value,
                          competing_demote = competing_demote),
            summary = list(
                n_supported = as.integer(result$call[["supported"]]),
                n_probable = as.integer(result$call[["probable"]]),
                n_undetermined = as.integer(result$call[["undetermined"]]),
                n_refuted = as.integer(result$call[["refuted"]]),
                independent_lines = length(unique_lines),
                independent_providers = independent,
                unavailable_providers = names(Filter(
                    function(p) identical(p$state, "unavailable"), provider_report))
            ),
            created_at = Sys.time()
        )
    )
    object <- sclet_log_command(
        object,
        "CallPhenotypeCells",
        params = list(name = name, providers = names(provider_report)),
        outputs = list(columns = paste0(prefix, c("call", "claim", "n_lines")))
    )

    for (pname in names(provider_report)) {
        pr <- provider_report[[pname]]
        if (identical(pr$state, "unavailable")) {
            warning(sprintf(
                "Evidence provider '%s' (%s) is unavailable, so it does not count toward the call: %s",
                pname, pr$label, paste(pr$issues, collapse = "; ")
            ), call. = FALSE)
        }
    }
    if (!is.null(competing_report)) {
        demoting <- names(Filter(function(x) isTRUE(x$demote), competing_report))
        if (length(demoting) > 0) {
            warning(sprintf(
                "Competing hypothesis '%s' separates the called cells from the rest beyond %.2f; calls that were supported or probable are demoted to undetermined because the split may track that programme rather than this phenotype.",
                demoting[[1]], competing_demote
            ), call. = FALSE)
        }
    }
    object
}

# ---------------------------------------------------------------------------
# registry plumbing
# ---------------------------------------------------------------------------

sclet_phenotype_registry_env <- new.env(parent = emptyenv())

sclet_phenotype_registry <- function() {
    if (!exists("registry", envir = sclet_phenotype_registry_env, inherits = FALSE)) {
        sclet_phenotype_registry_env$registry <- sclet_phenotype_builtin_providers()
    }
    get("registry", envir = sclet_phenotype_registry_env, inherits = FALSE)
}

sclet_phenotype_registry_set <- function(registry) {
    assign("registry", registry, envir = sclet_phenotype_registry_env)
    invisible(registry)
}

# The built-in providers. Each one is a thin adapter over an exported function so
# that its availability logic has exactly one definition.
sclet_phenotype_builtin_providers <- function() {
    list(
        cnv = list(
            name = "cnv",
            label = "Copy-number burden (inferCNV)",
            dependence = "independent",
            available = function(object, cnv_id = "cnv", ...) {
                rec <- sclet_get_analysis(object, "malignant_cnv")
                if (is.null(rec) || is.null(rec[[cnv_id]])) {
                    return(list(available = FALSE, issues = sprintf(
                        "No CNV record '%s'; run RunCNVInference().", cnv_id)))
                }
                rec <- rec[[cnv_id]]
                if (!identical(rec$status, "completed")) {
                    return(list(available = FALSE, issues = rec$reason %||%
                                    "CNV run did not complete"))
                }
                list(available = TRUE, issues = character(0))
            },
            decide = function(object, cnv_id = "cnv", cnv_cut = NULL,
                              reference_quantile = 0.99, ...) {
                cdata <- SummarizedExperiment::colData(object)
                burden <- cdata[[paste0(cnv_id, "_burden")]]
                is_ref <- cdata[[paste0(cnv_id, "_is_reference")]]
                ref_burden <- burden[!is.na(is_ref) & is_ref]
                cut <- if (!is.null(cnv_cut)) {
                    cnv_cut
                } else if (length(ref_burden) >= 5) {
                    sclet_cnv_reference_cut(ref_burden, reference_quantile)
                } else {
                    NA_real_
                }
                state <- if (is.na(cut)) {
                    rep("unavailable", nrow(cdata))
                } else {
                    ifelse(is.na(burden), "unavailable",
                           ifelse(burden > cut, "positive", "negative"))
                }
                list(state = state, detail = list(cut = cut, cnv_id = cnv_id))
            }
        ),
        module = list(
            name = "module",
            label = "Functional module score",
            # A module score is computed from the SAME expression matrix as the
            # signature, so it is not a second observation unless the caller
            # registers their own provider over an independent measurement.
            dependence = "shared_with_signature",
            available = function(object, module_genes = NULL, ...) {
                if (is.null(module_genes) || length(intersect(module_genes, rownames(object))) < 3) {
                    return(list(available = FALSE, issues =
                                    "Supply `module_genes` with at least 3 genes present in the object."))
                }
                list(available = TRUE, issues = character(0))
            },
            decide = function(object, module_genes = NULL, module_cut = NULL, ...) {
                genes <- intersect(module_genes, rownames(object))
                assay <- if ("logcounts" %in% SummarizedExperiment::assayNames(object)) "logcounts" else NULL
                if (is.null(assay)) {
                    return(list(state = NULL))
                }
                expr <- sclet_malignancy_dense(SummarizedExperiment::assay(object, assay))
                bins <- sclet_gene_expression_bins(expr, n_bins = 25)
                sc <- sclet_add_module_score(expr, genes, bins = bins, n_ctrl = 100,
                                             exclude = genes)
                cut <- if (is.null(module_cut)) stats::median(sc, na.rm = TRUE) else module_cut
                list(state = ifelse(is.na(sc), "unavailable",
                                    ifelse(sc >= cut, "positive", "negative")),
                     detail = list(cut = cut, genes = length(genes)))
            }
        )
    )
}

sclet_phenotype_resolve_signature <- function(signature) {
    if (is.null(signature)) {
        stop("Supply `signature` or `signature_score`.")
    }
    if (is.character(signature)) {
        return(signature)
    }
    if (is.list(signature) && !is.null(signature$tumor_genes)) {
        return(signature$tumor_genes)
    }
    if (is.list(signature) && !is.null(names(signature))) {
        idx <- match(c("malignant", "tumor"), tolower(names(signature)))
        idx <- idx[!is.na(idx)]
        if (length(idx) == 0) {
            stop("A named `signature` list needs a 'malignant' or 'tumor' element.")
        }
        return(signature[[idx[[1]]]])
    }
    stop("`signature` must be a gene vector, a named list, or FindMalignantSeedGenes() output.")
}
