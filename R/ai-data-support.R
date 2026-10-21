# M12-b: `analysis_story$data_support` -- what this data can and cannot answer.
#
# `analysis_story` records what the object has DONE (timeline, user decisions,
# conflicts, evidence gaps). It never said what the object can SUPPORT, so a
# planner had to infer that from assay names and column names -- which is
# exactly the kind of inference the 2026-10-16/19 sweeps measured going wrong
# (velocity planned on a counts-only object, milo planned with one sample per
# condition, a spatial question asked of a non-spatial object).
#
# Two things separate this from a second gate:
#
#   1. Every `reason` is an EXISTING gate name, taken from the readiness
#      diagnostics and the execution registry (`counts_assay_missing`,
#      `spliced_unspliced_missing`, `reference_missing`, `no_spatial_record`,
#      ...). A pin test asserts each reason below is a name the code already
#      uses, so this file cannot drift into a private vocabulary.
#   2. It never decides anything. It reports; ValidateAIPlan() and the action
#      prerequisites remain the authority. A `supported = TRUE` here does not
#      make an action runnable.
#
# The `unexplored` list is the part with no predecessor: `evidence_gaps` says
# "an analysis lacks evidence", `unexplored` says "the data was never analysed
# for this at all". That is the input for choosing what the next round should
# do, and it is computed from the object plus the ledger, never from a model.

# Gate names this file is allowed to emit. Kept as one vector so the pin test
# can check them against the code instead of against a comment.
sclet_ai_data_support_gates <- function() {
    c(
        "counts_assay_missing",
        "spliced_unspliced_missing",
        "reduction_missing",
        "reference_missing",
        "group_column_missing",
        "design_role_missing",
        "design_semantics_not_confirmed",
        "no_spatial_record",
        "no_multimodal_record",
        "cluster_identity_missing",
        "terminal_state_codes_missing",
        "velocity_record_missing"
    )
}

# Readiness questions are phrased for a person; `blocked` needs the gate name.
# Mapping the action to the gate that governs it keeps `data_support` and the
# action prerequisites saying the same thing in the same words.
sclet_ai_data_support_gate_for <- function(action, checks) {
    if (identical(action, "run_velocity")) {
        if (isFALSE(checks$spliced_unspliced_available)) return("spliced_unspliced_missing")
        if (isFALSE(checks$counts_available)) return("counts_assay_missing")
        return(NULL)
    }
    if (identical(action, "run_fate_analysis")) {
        if (isFALSE(checks$velocity_record_available)) return("velocity_record_missing")
        return(NULL)
    }
    if (identical(action, "run_milo")) {
        if (isFALSE(checks$design_roles_confirmed)) return("design_semantics_not_confirmed")
        if (isFALSE(checks$sample_column_resolved) || isFALSE(checks$condition_column_resolved)) {
            return("design_role_missing")
        }
        return(NULL)
    }
    if (identical(action, "run_annotation")) {
        if (isFALSE(checks$reference_available)) return("reference_missing")
        return(NULL)
    }
    if (identical(action, "run_cci")) {
        if (isFALSE(checks$method_available)) return("reference_missing")
        return(NULL)
    }
    NULL
}

# The shapes block: facts about the data, each one a boolean or a list of
# candidate column names. Nothing here is a verdict.
sclet_ai_data_support_shapes <- function(object, columns, assay_names) {
    has_assay <- function(name) name %in% assay_names
    n_cells <- ncol(object)

    # Counts must be checked for integrality, not just for the name "counts":
    # a log-normalized matrix stored under `counts` is still not counts, and
    # velocity/inferCNV need the real thing.
    counts_integer <- tryCatch({
        if (!has_assay("counts")) {
            FALSE
        } else {
            values <- SummarizedExperiment::assay(object, "counts")
            if (inherits(values, "dgCMatrix") || inherits(values, "dgTMatrix") ||
                inherits(values, "matrix")) {
                sample_values <- if (inherits(values, "matrix")) {
                    values[seq_len(min(50L, nrow(values))), seq_len(min(50L, ncol(values))), drop = FALSE]
                } else {
                    as.numeric(utils::head(as.vector(values@x), 200L))
                }
                sample_values <- sample_values[!is.na(sample_values)]
                length(sample_values) > 0L &&
                    all(sample_values >= 0) &&
                    all(abs(sample_values - round(sample_values)) < 1e-8)
            } else {
                TRUE
            }
        }
    }, error = function(e) FALSE)

    # Column shapes are read by NAME PATTERN and then reported as counts only.
    # No level labels leave this function: a donor id or a patient code is
    # identifying, a number of levels is not.
    level_count <- function(column) {
        tryCatch({
            values <- as.character(SummarizedExperiment::colData(object)[[column]])
            length(unique(values[!is.na(values)]))
        }, error = function(e) NA_integer_)
    }
    pattern_columns <- function(pattern) {
        columns[grepl(pattern, columns, ignore.case = TRUE)]
    }
    sample_columns <- union(
        pattern_columns("sample"),
        pattern_columns("donor|patient|subject|individual")
    )
    condition_columns <- pattern_columns("condition|treatment|stim|group|genotype|time")
    batch_columns <- pattern_columns("batch|lane|run|plate|pool")

    multi_sample <- any(vapply(sample_columns, function(column) {
        n <- level_count(column)
        !is.na(n) && n > 1L
    }, logical(1L)))

    # Replication means: at least one condition has more than one sample. It is
    # computed from the pair, so it stays FALSE when either side is unknown --
    # "unknown" must not read as "no".
    #
    # Written with a plain assignment inside the tryCatch expression, not with
    # `return(TRUE)`: a `return()` there returns from the ENCLOSING function
    # (`sclet_ai_data_support_shapes`), so the first object with replicable
    # samples made the whole shapes block evaluate to a bare TRUE. Caught by
    # the determinism pin, which is why that pin exists.
    has_replicates <- FALSE
    if (length(sample_columns) && length(condition_columns)) {
        has_replicates <- tryCatch({
            cols <- SummarizedExperiment::colData(object)
            found <- FALSE
            for (s in sample_columns) {
                for (cc in condition_columns) {
                    pairing <- unique(data.frame(
                        sample = as.character(cols[[s]]),
                        condition = as.character(cols[[cc]]),
                        stringsAsFactors = FALSE
                    ))
                    pairing <- pairing[!is.na(pairing$sample) & !is.na(pairing$condition), , drop = FALSE]
                    if (nrow(pairing) && any(table(pairing$condition) > 1L)) {
                        found <- TRUE
                        break
                    }
                }
                if (found) break
            }
            found
        }, error = function(e) FALSE)
    }

    # Reuse the package's own spacious-coordinate probe rather than calling
    # SpatialExperiment::spatialCoords directly: that accessor is not exported
    # from SingleCellExperiment, and R CMD check reports the "::" call as a
    # missing object unless it sits behind requireNamespace(). The existing
    # helper already handles the class check, the optional dependency and the
    # fallbacks (a `spatial` reduction, coordinate-like colData names).
    has_spatial <- tryCatch(
        sclet_ai_spatial_coordinates(object),
        error = function(e) FALSE
    )
    alternative_experiments <- tryCatch(
        SingleCellExperiment::altExpNames(object),
        error = function(e) character()
    )

    list(
        n_cells = as.integer(n_cells),
        n_genes = as.integer(nrow(object)),
        # M13-a: whether this object carries ambient-RNA pressure. Deterministic,
        # bounded, and explicitly a screen rather than a measurement -- see
        # sclet_ai_ambient_pressure().
        ambient = sclet_ai_ambient_pressure(object),
        assays = assay_names,
        has_counts = has_assay("counts"),
        counts_are_integer = counts_integer,
        has_spliced_unspliced = all(c("spliced", "unspliced") %in% assay_names),
        has_multiple_samples = multi_sample,
        has_replicates = has_replicates,
        has_spatial = has_spatial,
        has_multimodal = length(alternative_experiments) > 0L,
        candidate_sample_columns = sample_columns,
        candidate_condition_columns = condition_columns,
        candidate_batch_columns = batch_columns,
        candidate_design_columns = sort(unique(c(
            sample_columns, condition_columns, batch_columns
        ))),
        columns = columns,
        raw_values_included = FALSE
    )
}

# `supports`: per question class, whether the data can answer it, why not, and
# what would change that. Each `reason` comes from sclet_ai_data_support_gates().
sclet_ai_data_support_supports <- function(object, shapes, records = list()) {
    entry <- function(question_class, supported, reason = NULL, unblocked_by = NULL) {
        list(
            question_class = question_class,
            supported = isTRUE(supported),
            reason = reason,
            unblocked_by = unblocked_by
        )
    }

    descriptor_present <- function(action) {
        registry <- tryCatch(AIExecutionRegistry(object), error = function(e) NULL)
        !is.null(registry) && !is.null(registry[[action]])
    }

    supports <- list()

    # Velocity / fate: needs spliced+unspliced, and the counts assay that
    # velociraptor reads. Measured 2026-10-19 (M9): every velocity export in
    # this package carries only logcounts, so this gate is what stopped
    # planning rather than letting the run die inside Python.
    if (shapes$has_spliced_unspliced && shapes$counts_are_integer) {
        supports[[length(supports) + 1L]] <- entry("velocity", TRUE)
    } else if (!shapes$has_spliced_unspliced) {
        supports[[length(supports) + 1L]] <- entry(
            "velocity", FALSE, "spliced_unspliced_missing",
            "use a velocity export that carries spliced and unspliced assays, or omit velocity from the plan"
        )
    } else {
        supports[[length(supports) + 1L]] <- entry(
            "velocity", FALSE, "counts_assay_missing",
            "velociraptor reads counts; provide the unspliced/counts export rather than logcounts only"
        )
    }

    # Trajectory: needs an embedding and an identity, both of which the
    # standard spine produces. It does not need counts.
    #
    # "has a grouping" is read from the OBJECT and from any analysis that
    # presupposes one, NOT from a record literally typed "clustering". Measured
    # on the M11 pancreas object: it carries preprocess + reduction + velocity +
    # trajectory + cellrank and NO `clustering` record, because the runner went
    # straight from the embedding to the trajectory. Keying on the record type
    # reported `cluster_identity_missing` on an object that had a trajectory --
    # a false refusal, and it also silenced every `unexplored` entry that
    # required clustering, so the whole block came back empty.
    reductions <- tryCatch(SingleCellExperiment::reducedDimNames(object),
        error = function(e) character())
    grouping_evidence <- sclet_ai_data_support_grouping_evidence(object, records)
    if (length(reductions) && grouping_evidence$available) {
        supports[[length(supports) + 1L]] <- entry("trajectory", TRUE)
    } else if (!length(reductions)) {
        supports[[length(supports) + 1L]] <- entry(
            "trajectory", FALSE, "reduction_missing",
            "run the dimension-reduction steps first (run_pca, find_neighbors)"
        )
    } else {
        supports[[length(supports) + 1L]] <- entry(
            "trajectory", FALSE, "cluster_identity_missing",
            "cluster the cells first; a trajectory needs a grouping to order"
        )
    }

    # Differential abundance: the replication signal, not the package, is the
    # usual gate. `has_replicates` is FALSE both when there are no replicates
    # and when the design columns are unknown; the reason distinguishes them.
    if (shapes$has_replicates) {
        supports[[length(supports) + 1L]] <- entry("differential_abundance", TRUE)
    } else if (!length(shapes$candidate_sample_columns) ||
        !length(shapes$candidate_condition_columns)) {
        supports[[length(supports) + 1L]] <- entry(
            "differential_abundance", FALSE, "design_role_missing",
            "declare which columns are the sample and the condition (ConfirmAIDesignSemantics)"
        )
    } else {
        supports[[length(supports) + 1L]] <- entry(
            "differential_abundance", FALSE, "design_role_missing",
            "differential abundance attributes neighbourhood counts to samples; it needs at least two samples per condition"
        )
    }

    # Spatial and multimodal are read-only refusals by design.
    if (isTRUE(shapes$has_spatial)) {
        supports[[length(supports) + 1L]] <- entry("spatial", TRUE)
    } else {
        supports[[length(supports) + 1L]] <- entry(
            "spatial", FALSE, "no_spatial_record",
            "a spatial question needs spatial coordinates on the object"
        )
    }
    if (isTRUE(shapes$has_multimodal)) {
        supports[[length(supports) + 1L]] <- entry("multimodal", TRUE)
    } else {
        supports[[length(supports) + 1L]] <- entry(
            "multimodal", FALSE, "no_multimodal_record",
            "a multimodal question needs a second modality in altExp"
        )
    }

    # Ambient-RNA cleaning (M13-a). Available when there are counts to estimate
    # the profile from and the celda backend is installed. When the screen says
    # the data carries ambient pressure, this is the route that acts on it.
    has_integer_counts <- isTRUE(shapes$counts_are_integer)
    if (has_integer_counts && base::requireNamespace("celda", quietly = TRUE)) {
        supports[[length(supports) + 1L]] <- entry("ambient_cleaning", TRUE)
    } else if (!has_integer_counts) {
        supports[[length(supports) + 1L]] <- entry(
            "ambient_cleaning", FALSE, "counts_assay_missing",
            "decontX estimates the ambient profile from raw integer counts; provide a counts assay"
        )
    } else {
        supports[[length(supports) + 1L]] <- entry(
            "ambient_cleaning", FALSE, "optional_package_missing",
            "install celda to run decontX"
        )
    }

    # The descriptive spine is available on any object with an assay.
    supports[[length(supports) + 1L]] <- entry(
        "descriptive_spine",
        length(shapes$assays) > 0L,
        if (length(shapes$assays)) NULL else "reduction_missing"
    )

    supports
}

# `unexplored`: data value the ledger has no result for. This is the list that
# answers "what should the next round do", and it is the field with no
# predecessor in `analysis_story`.
#
# The rule for including an entry is deliberately conservative: only value that
# the data can objectively support AND that no record covers. A suggestion that
# the data cannot support would be worse than silence, because a planner trusts
# it.
# Does the object have a cell grouping? Read from every source that can prove
# one, because a record named "clustering" is only the most obvious of them.
#
# Measured 2026-10-21 on the M11 pancreas object: it has preprocess, reduction,
# velocity, trajectory and cellrank records and no `clustering` record, because
# that run went from the embedding straight to a trajectory. Keying grouping on
# the record type alone made two wrong claims at once --
# `cluster_identity_missing` on an object that already had a trajectory, and an
# empty `unexplored` list, because every suggestion there required clustering.
sclet_ai_data_support_grouping_evidence <- function(object, records = list()) {
    analyses <- if (is.list(records$analyses)) records$analyses else list()
    types <- vapply(analyses, function(record) {
        if (is.list(record)) as.character(record$type %||% "") else ""
    }, character(1))

    # An analysis that cannot be planned without a grouping is itself proof
    # that one existed when it ran.
    presupposing <- c("clustering", "detest", "enrichment", "annotation",
        "trajectory", "velocity", "cellrank", "milo", "rare_cells",
        "doublet_detection", "cnv", "malignant_cnv", "communication", "cci")

    sources <- character()
    if ("clustering" %in% types) sources <- c(sources, "clustering_record")
    implied <- intersect(presupposing, types)
    implied <- setdiff(implied, "clustering")
    if (length(implied)) {
        sources <- c(sources, paste0("implied_by:", paste(sort(implied), collapse = "+")))
    }
    # The object itself can carry a grouping without any AI record.
    has_ident <- tryCatch(!is.null(object$sclet_ident), error = function(e) FALSE)
    if (isTRUE(has_ident)) sources <- c(sources, "object_sclet_ident")
    has_labels <- tryCatch({
        length(unique(as.character(object$label[!is.na(object$label)]))) > 1L
    }, error = function(e) FALSE)
    if (isTRUE(has_labels)) sources <- c(sources, "object_labels")

    list(
        available = length(sources) > 0L,
        sources = sources,
        record_types = sort(unique(types[nzchar(types)]))
    )
}

sclet_ai_data_support_unexplored <- function(object, shapes, supports, analyses, records = list()) {
    has_analysis_type <- function(type) {
        if (!is.list(analyses)) return(FALSE)
        any(vapply(analyses, function(record) {
            is.list(record) && identical(as.character(record$type %||% ""), type)
        }, logical(1L)))
    }
    supported <- function(class) {
        any(vapply(supports, function(x) {
            identical(x$question_class, class) && isTRUE(x$supported)
        }, logical(1L)))
    }
    # A grouping is proved by any analysis that presupposes one, not by a
    # record literally typed "clustering" -- see
    # sclet_ai_data_support_grouping_evidence() for the measured case.
    grouping <- sclet_ai_data_support_grouping_evidence(object, records)
    has_grouping <- function() isTRUE(grouping$available)

    entries <- list()
    add <- function(what, action, why) {
        entries[[length(entries) + 1L]] <<- list(
            what = what,
            action = action,
            why = why
        )
    }

    # A key omission on purpose: `unexplored` NEVER suggests re-running an
    # analysis the object already has. "Tested and found nothing" and "never
    # tested" are different states, and only the second belongs here.

    if (isTRUE(shapes$has_multiple_samples) &&
        length(shapes$candidate_batch_columns) &&
        !has_analysis_type("integration") &&
        !has_analysis_type("integration_routes")) {
        add(
            "colData names sample and batch columns but no integration has been run, so batch structure has not been assessed",
            "run_integration",
            "several samples with a batch column is the shape integration exists for"
        )
    }

    if (supported("differential_abundance") &&
        !has_analysis_type("milo") &&
        has_grouping()) {
        add(
            "the design has replicate samples per condition and a clustering exists, but no differential abundance has been computed",
            "run_milo",
            "replication is what makes neighbourhood abundance testable"
        )
    }

    if (supported("velocity") && !has_analysis_type("velocity")) {
        add(
            "the object carries spliced/unspliced and integer counts, so velocity is available and has not been run",
            "run_velocity",
            "the assays velocity needs are present"
        )
    }

    if (supported("trajectory") && !has_analysis_type("trajectory")) {
        add(
            "an embedding and a grouping exist but no trajectory has been ordered",
            "run_trajectory",
            "ordering clusters along a lineage is available on this object"
        )
    }

    if (has_grouping() && !has_analysis_type("detest") &&
        !has_analysis_type("markers")) {
        add(
            "cells are clustered but no marker test has been run, so the clusters have no supporting genes",
            "run_de_test",
            "a cluster without markers is a grouping, not a finding"
        )
    }

    # Ambient first, because it changes how every downstream number should be
    # read: suggesting markers on contaminated data invites interpreting ambient
    # signal as biology.
    # RunDecontX() writes a `preprocess` record carrying id "decontX". Checking
    # the TYPE alone would be wrong: every normalisation is a preprocess record
    # too, so "preprocess exists" says nothing about whether cleaning happened.
    has_decontx <- FALSE
    if (is.list(analyses)) {
        for (record in analyses) {
            if (!is.list(record)) next
            if (!identical(as.character(record$type %||% ""), "preprocess")) next
            record_id <- tolower(as.character(record$id %||% ""))
            if (grepl("decontx", record_id, fixed = TRUE)) has_decontx <- TRUE
        }
    }
    ambient <- shapes$ambient %||% list()
    if (isTRUE(ambient$available) && identical(ambient$label, "high") &&
        supported("ambient_cleaning") && !has_decontx) {
        add(
            paste0(
                "the ambient-pressure screen is high (low-UMI tail share ",
                ambient$low_umi_tail_share, ", top-gene detection uniformity ",
                ambient$top_gene_detection_uniformity,
                "), and no ambient-RNA cleaning has been run"
            ),
            "run_decontx",
            "a high screen means contamination is worth removing before the composition and the marker tables are interpreted"
        )
    }

    if (has_analysis_type("detest") && !has_analysis_type("enrichment")) {
        add(
            "a marker table exists but no enrichment was run on it",
            "run_enrichment",
            "enrichment interprets the marker table that already exists"
        )
    }

    if (has_grouping() && !has_analysis_type("annotation")) {
        add(
            "cells are clustered but no reference annotation has been applied",
            "run_annotation",
            "label transfer is available once clusters exist; the reference is the user's choice"
        )
    }

    if (supported("spatial") && !has_analysis_type("spatial_deconvolution")) {
        add(
            "the object carries spatial coordinates but no spatial analysis has been run",
            "run_spatial_deconvolution",
            "spatial coordinates are present"
        )
    }

    entries
}

# Build the whole block. `records` is what sclet_ai_build_analysis_story()
# already holds, so this adds no second pass over the state store.
sclet_ai_build_data_support <- function(object, records) {
    analyses <- if (is.list(records$analyses)) records$analyses else list()
    columns <- tryCatch(
        colnames(SummarizedExperiment::colData(object)),
        error = function(e) character()
    )
    assay_names <- tryCatch(
        SummarizedExperiment::assayNames(object),
        error = function(e) character()
    )

    shapes <- sclet_ai_data_support_shapes(object, columns, assay_names)
    supports <- sclet_ai_data_support_supports(object, shapes, records)

    blocked <- list()
    for (s in supports) {
        if (isTRUE(s$supported)) next
        blocked[[length(blocked) + 1L]] <- list(
            question_class = s$question_class,
            gate = s$reason,
            unblocked_by = s$unblocked_by
        )
    }

    # Unconfirmed semantics: the choices the planner must not default. Only
    # choices the data MAKES relevant are listed -- an annotation reference is
    # not "unconfirmed" on an object that was never clustered.
    unconfirmed <- list()
    if (length(shapes$candidate_sample_columns) && !length(shapes$candidate_condition_columns)) {
        unconfirmed[[length(unconfirmed) + 1L]] <- list(
            choice = "condition",
            reason = "a sample-like column exists but no condition-like column does",
            gate = "design_role_missing"
        )
    }
    if (length(shapes$candidate_design_columns)) {
        unconfirmed[[length(unconfirmed) + 1L]] <- list(
            choice = "batch_and_condition_roles",
            reason = "column names suggest design roles but no ConfirmAIDesignSemantics() declaration was made",
            gate = "design_semantics_not_confirmed"
        )
    }
    if (length(shapes$assays) && any(vapply(analyses, function(record) {
        is.list(record) && identical(as.character(record$type %||% ""), "clustering")
    }, logical(1L))) && !any(vapply(analyses, function(record) {
        is.list(record) && identical(as.character(record$type %||% ""), "annotation")
    }, logical(1L)))) {
        unconfirmed[[length(unconfirmed) + 1L]] <- list(
            choice = "annotation_reference",
            reason = "clusters exist and no annotation has been applied; the reference is the user's choice, never a default",
            gate = "reference_missing"
        )
    }

    list(
        schema_version = "1.0",
        shapes = shapes,
        supports = supports,
        blocked = blocked,
        unconfirmed_semantics = unconfirmed,
        unexplored = sclet_ai_data_support_unexplored(object, shapes, supports, analyses, records),
        raw_values_included = FALSE
    )
}

# ---- M13-b: does this round respond to the previous one? --------------------
#
# D3 of the four plan-quality dimensions, made computable. `data_support` (§ above)
# answers "what can this data do"; this answers "did the planner ACT on what the
# last round said".
#
# The input is two rounds of the same object: the previous plan and report, and the
# current plan. Nothing here reads the model's reasoning -- that would be
# unfalsifiable. Every judge below is a relation between the two plans and the
# report's own fields.
#
# Four judges, and the one that has to exist: `regressed`. A dimension that can only
# report "good" or "average" cannot report "worse", and "this round undid what the
# last round established" is the single most expensive progress a plan can make.

# The actions that remove a warning, by the kind of warning they answer. Only
# pairs where the fix is an actual action are listed -- a warning with no
# answering action stays "carried", which is a reportable state, not a gap.
sclet_ai_round_progression_remedies <- function() {
    # EVERY action named here must be a registry key, and a pin enforces it.
    #
    # History worth keeping: this table first named `run_decontx` and
    # `run_imputation` from memory. Neither existed -- RunDecontX() was exported
    # but had no AIAction, and there is no imputation action at all. A remedy
    # naming a non-existent action never matches, so an ambient warning read as
    # "ignored" forever, and the judge was silently wrong with no symptom. A
    # real-data probe returned `regressed` and the reason had to be traced by
    # hand. M13-a then created the `run_decontx` action, so `ambient` has a
    # remedy for the first time. There is still no imputation action, and this
    # table does not pretend otherwise.
    list(
        doublet = c("run_doublet_detection"),
        ambient = c("run_decontx"),
        contamination = c("run_decontx"),
        batch = c("run_integration", "confirm_design"),
        design = c("confirm_design"),
        reference = c("run_annotation"),
        annotation = c("run_annotation"),
        coverage = c("run_enrichment", "run_de_test"),
        marker = c("run_de_test"),
        enrichment = c("run_enrichment")
    )
}

# The subject of a warning, as a token that can be matched against a remedy.
# A warning is a sentence; this reduces it to what it is ABOUT, which is the only
# part an action can respond to. Free text the whole way through would make the
# judge a keyword lottery.
sclet_ai_round_progression_warning_topic <- function(warning) {
    w <- tolower(warning)
    if (!nzchar(trimws(w))) return(NULL)
    for (topic in c("doublet", "ambient", "contaminat", "batch", "design",
        "reference", "annotation", "coverage", "marker", "enrichment")) {
        if (grepl(topic, w, fixed = TRUE)) return(topic)
    }
    NULL
}

sclet_ai_round_progression <- function(
    previous_plan,
    previous_report,
    current_plan,
    current_report = NULL
) {
    # ---- normalise the two plans into comparable shape -------------------
    actions_of <- function(plan) {
        if (is.null(plan$actions) || !length(plan$actions)) {
            return(character())
        }
        vapply(plan$actions, function(x) as.character(x$action %||% ""), character(1))
    }
    params_of <- function(plan) {
        out <- list()
        for (step in (plan$actions %||% list())) {
            key <- paste0(as.character(step$action %||% ""), "#", as.character(step$id %||% ""))
            out[[key]] <- step$params %||% list()
        }
        out
    }
    prev_actions <- actions_of(previous_plan)
    curr_actions <- actions_of(current_plan)
    prev_params <- params_of(previous_plan)
    curr_params <- params_of(current_plan)

    # ---- D3.2 no-progress: same action, same params, repeated ------------
    # A step is "repeated" only when the action AND its parameters are unchanged.
    # The same action with new parameters is a new analysis, not a replay.
    prev_by_action <- split(names(prev_params), vapply(names(prev_params),
        function(k) sub("#.*$", "", k), character(1)))
    repeated <- character()
    for (action in intersect(names(prev_by_action), unique(curr_actions))) {
        old_keys <- prev_by_action[[action]]
        new_keys <- curr_actions %in% action
        new_names <- names(curr_params)[new_keys]
        for (nk in new_names) {
            # find a previous step of the same action with equal params
            for (ok in old_keys) {
                if (identical(prev_params[[ok]], curr_params[[nk]])) {
                    repeated <- c(repeated, nk)
                    break
                }
            }
        }
    }
    progress_ratio <- if (length(curr_actions)) {
        1 - length(repeated) / length(curr_actions)
    } else {
        NA_real_
    }

    # ---- D3.1 answered warnings -----------------------------------------
    # Each warning reduces to a topic; a topic is answered if this round runs an
    # action that remediates it, or if the rationale states it is carried forward
    # deliberately. "Deliberately carried" is an answer, not silence: a warning the
    # user has accepted and recorded is not a defect.
    warnings <- character()
    raw <- previous_report$warnings %||% list()
    if (is.list(raw)) {
        warnings <- vapply(raw, function(x) as.character(x)[1L], character(1))
    } else if (is.character(raw)) {
        warnings <- as.character(raw)
    }
    warnings <- warnings[nzchar(warnings)]
    rationale <- if (is.null(current_plan$rationale)) {
        ""
    } else {
        paste(as.character(current_plan$rationale), collapse = " ")
    }
    rationale_l <- tolower(rationale)

    remedies <- sclet_ai_round_progression_remedies()
    n_answered <- 0L
    n_carried <- 0L
    n_ignored <- 0L
    unanswered <- character()
    for (warning in warnings) {
        topic <- sclet_ai_round_progression_warning_topic(warning)
        answered <- FALSE
        carried <- FALSE
        if (!is.null(topic)) {
            candidates <- remedies[[topic]] %||% character()
            if (any(candidates %in% curr_actions)) answered <- TRUE
        }
        # A carried warning must be NAMED, not merely unmentioned. "We cannot fix
        # X because Y" counts; saying nothing about X does not.
        # `nzchar(NULL)` is a runtime error, and a warning with no recognisable
        # topic is exactly the kind of input that produces one -- so test the
        # shape before the call, not inside it.
        if (!answered && length(topic) == 1L && !is.na(topic) && nzchar(topic) &&
            grepl(topic, rationale_l, fixed = TRUE)) {
            carried <- TRUE
        }
        if (answered) {
            n_answered <- n_answered + 1L
        } else if (carried) {
            n_carried <- n_carried + 1L
        } else {
            n_ignored <- n_ignored + 1L
            unanswered <- c(unanswered, warning)
        }
    }

    # ---- D3.3 evidence gaps closed --------------------------------------
    gaps <- previous_report$evidence_gaps %||% list()
    n_gaps <- length(gaps)
    n_gaps_closed <- 0L
    if (n_gaps) {
        for (gap in gaps) {
            if (!is.list(gap)) next
            analysis_type <- as.character(gap$analysis_type %||% "")
            reason <- as.character(gap$reason %||% gap$gap %||% "")
            closed <- FALSE
            # A gap names an analysis type (rare_cells) and often the signal that
            # is missing (doublet_evidence_available): either this round produces
            # it, or the analysis itself now runs.
            if (nzchar(analysis_type) && any(grepl(analysis_type, curr_actions))) {
                closed <- TRUE
            }
            if (!closed && nzchar(reason)) {
                token <- sub("_.*$", "", reason)
                if (nzchar(token) && any(grepl(token, curr_actions))) closed <- TRUE
            }
            if (closed) n_gaps_closed <- n_gaps_closed + 1L
        }
    }

    # ---- D3.4 recorded user decisions actually used ----------------------
    decisions <- previous_report$user_decisions %||% list()
    n_decisions <- length(decisions)
    n_decisions_used <- 0L
    unused_decisions <- character()
    for (decision in decisions) {
        if (!is.list(decision)) next
        summary <- as.character(decision$summary %||% "")
        if (!nzchar(summary)) next
        # A decision is used if this round names its subject. The clarification
        # channel records the answer as a short phrase ("use MonacoImmuneData as
        # the annotation reference"), and matching its FIRST token matches the
        # verb, not the subject -- so "use" never appears in the next plan and
        # every recorded decision reads as unused. A stop-word list is applied
        # and the longest remaining word is the subject: "monacoimmunedata".
        # The split class must keep case while matching case-insensitively, or
        # `MonacoImmuneData` is cut at every capital and the "subject word" comes
        # out as "onaco" / "mmune" / "ata". Split on characters that are not
        # letters or digits in ANY case.
        words <- strsplit(summary, "[^A-Za-z0-9._]+")[[1L]]
        words <- words[nzchar(words)]
        words <- tolower(words)
        words <- setdiff(words, c(
            "use", "using", "as", "the", "a", "an", "to", "for", "of", "with",
            "chosen", "choice", "user", "decision", "recorded", "selected",
            "is", "be", "by", "in", "on", "and", "or", "not"
        ))
        token <- if (length(words)) {
            words[which.max(nchar(words))]
        } else {
            ""
        }
        used <- nzchar(token) && nchar(token) >= 6L && (
            grepl(token, rationale_l, fixed = TRUE) ||
            any(grepl(token, curr_actions, ignore.case = TRUE))
        )
        if (used) {
            n_decisions_used <- n_decisions_used + 1L
        } else {
            unused_decisions <- c(unused_decisions, summary)
        }
    }

    # ---- regressed: this round undoes an established answer --------------
    # The case with no positive counterpart elsewhere in this block. A round that
    # re-analyses a population whose contamination the previous round already
    # established is worse than a round that does nothing, because it spends
    # effort re-deriving a discarded answer.
    established_contamination <- FALSE
    # vapply needs a length-1 result; the topic helper returns NULL for a warning
    # with no recognisable subject, which is length 0 and aborts the whole call.
    # Wrapping makes "no topic" an empty string instead of a crash.
    prev_topics <- vapply(warnings, function(w) {
        t <- sclet_ai_round_progression_warning_topic(w)
        if (length(t) == 1L && !is.na(t)) t else ""
    }, character(1))
    prev_topics <- prev_topics[nzchar(prev_topics)]
    for (topic in unique(prev_topics)) {
        # the previous round warned about X, and this round runs an action whose
        # output is exactly the thing being warned about
        reanalysing <- c(
            identical(topic, "doublet") && any(c(
                "run_annotation", "run_de_test", "run_cci"
            ) %in% curr_actions) && !any(c(
                "run_doublet_detection", "run_decontx"
            ) %in% curr_actions)
        )
        if (any(reanalysing)) established_contamination <- TRUE
    }

    # ---- the class -------------------------------------------------------
    progression_class <- "advanced"
    if (established_contamination) {
        progression_class <- "regressed"
    } else if (!is.na(progress_ratio) && progress_ratio == 0 &&
        !n_answered && !n_gaps_closed && !n_decisions_used) {
        progression_class <- "no_progress"
    } else if (!n_answered && !n_gaps_closed && !n_decisions_used) {
        progression_class <- "confirmatory"
    }

    list(
        schema_version = "1.0",
        # the four judges, flat so a record can cite them directly
        n_warnings = length(warnings),
        n_warnings_answered = n_answered,
        n_warnings_carried = n_carried,
        n_warnings_ignored = n_ignored,
        unanswered_warnings = unanswered,
        n_actions = length(curr_actions),
        n_repeated_steps = length(repeated),
        repeated_steps = repeated,
        progress_ratio = progress_ratio,
        n_evidence_gaps = n_gaps,
        n_evidence_gaps_closed = n_gaps_closed,
        n_user_decisions = n_decisions,
        n_user_decisions_used = n_decisions_used,
        unused_user_decisions = unused_decisions,
        regressed = isTRUE(established_contamination),
        progression_class = progression_class,
        raw_values_included = FALSE
    )
}

# ---- M13-b: find the previous round's plan and report on this object -------
#
# `sclet_ai_round_progression()` judges a pair of rounds; this finds the pair.
#
# The previous round is identified from the ledger, not from a parameter, because a
# caller that has to remember to pass it will eventually not -- and then D3 silently
# becomes "always a first round", which is indistinguishable from a workflow that
# never advances. Locating it from the ledger makes the failure visible as
# `no previous round` instead.
sclet_ai_previous_round_progression <- function(object, current_plan, current_report) {
    no_previous <- list(
        schema_version = "1.0",
        reason = "no previous round: this is the first round on this object",
        progression_class = "first_round",
        raw_values_included = FALSE
    )
    # NOTE: no `return()` anywhere inside the tryCatch expression. A `return()`
    # there returns from THIS function (sclet_ai_previous_round_progression),
    # not from the tryCatch -- measured 2026-10-21, when an object with no
    # recorded report made this helper hand back a bare NULL and every D3 field
    # vanished. Every early exit assigns to a variable and breaks out with
    # `break`/a value instead. This is the third instance of the same shape in
    # this package; the pin below exists because the shape keeps coming back.
    previous <- NULL
    # The previous round is assembled from TWO records, because no single record
    # holds both halves:
    #
    #   - the plan (what was run) comes from the newest `ai_execution` record's
    #     `params$steps`, which lists every action with its resolved params;
    #   - the report (what the run then said) comes from the newest `ai_report`
    #     record, whose `artifacts$warnings` is the warnings list and whose
    #     `summary$cited_user_decisions` is what the report itself cited.
    #
    # Verified against the real M1-M11 objects: reports do NOT carry a plan, and
    # `warnings` sit in artifacts, not summary. Looking for either in the other
    # place produced `not_comparable` on every real object -- a silent failure
    # that reads like "D3 was never computed here", which is indistinguishable
    # from a workflow that never advances.
    found <- tryCatch({
        ledger <- GetAnalysisLedger(object, detail = "full", include_artifacts = TRUE)
        keys <- names(ledger$analyses %||% list())
        is_execution <- vapply(keys, function(k) {
            identical(as.character(ledger$analyses[[k]]$type %||% ""), "ai_execution")
        }, logical(1L))
        is_report <- vapply(keys, function(k) {
            grepl("^ai_report", k) || identical(
                as.character(ledger$analyses[[k]]$type %||% ""), "ai_ai_report"
            )
        }, logical(1L))
        exec_keys <- keys[is_execution]
        report_keys <- keys[is_report]
        if (!length(exec_keys) || !length(report_keys)) {
            character(0)
        } else {
            execution <- ledger$analyses[[exec_keys[[length(exec_keys)]]]]
            report <- ledger$analyses[[report_keys[[length(report_keys)]]]]
            list(execution = execution, report = report)
        }
    }, error = function(e) character(0))
    if (!identical(found, character(0)) && !is.null(found)) {
        previous <- found
    }
    if (is.null(previous)) return(no_previous)

    # Field locations on a recorded AI report, verified against the real M1-M11
    # objects rather than assumed: `warnings` live in ARTIFACTS (as a list of
    # character vectors), and the user decisions the report cited are
    # `summary$cited_user_decisions`. Reading them off `summary$warnings` --
    # where they are not -- made every round look silencingly empty.
    execution <- previous$execution %||% list()
    report <- previous$report %||% list()
    summary <- report$summary %||% list()
    artifacts <- report$artifacts %||% list()

    # Rebuild the previous round's plan from its execution record's steps, so the
    # judges compare two plans rather than a plan against nothing.
    steps <- (execution$params %||% list())$steps %||% list()
    prev_plan <- if (length(steps)) {
        new_sclet_ai_plan(
            task = paste0("m13b_prev_", sub("^ai_execution_plan_", "",
                as.character(execution$id %||% "unknown"))),
            actions = lapply(steps, function(step) {
                list(
                    id = as.character(step$id %||% ""),
                    action = as.character(step$action %||% ""),
                    params = step$params %||% list()
                )
            }),
            context_fingerprint = as.character(
                (execution$inputs %||% list())$context_fingerprint %||% NA_character_
            )
        )
    } else {
        NULL
    }
    prev_report <- list(
        warnings = artifacts$warnings %||% summary$warnings %||% list(),
        evidence_gaps = summary$evidence_gaps %||% list(),
        user_decisions = summary$cited_user_decisions %||% list()
    )
    if (is.null(prev_plan)) {
        # a report with no plan behind it cannot be compared: say so rather than
        # comparing against an empty plan, which would make every round look like
        # progress
        return(list(
            schema_version = "1.0",
            reason = "the previous round recorded a report but no plan to compare against",
            progression_class = "not_comparable",
            raw_values_included = FALSE
        ))
    }

    sclet_ai_round_progression(
        prev_plan,
        prev_report,
        current_plan,
        current_report
    )
}

# ---- M13-a: ambient pressure, computed from the object ---------------------
#
# A deterministic scalar saying whether this object carries ambient-RNA pressure,
# so that "the data is dirty in this specific way" is a fact in the context rather
# than something a planner has to infer from a contamination-looking cluster.
#
# TWO heuristics, reported SEPARATELY, and deliberately not multiplied into one
# number. Multiplying two soft signals produces a figure with two decimal places
# and no meaning; a reader cannot tell which half moved. Each is 0..1 and each
# says what it measures.
#
# The thresholds below are CALIBRATED ON THE FIXTURES THIS PACKAGE RUNS, not
# derived from the literature. They are reported alongside the numbers so a reader
# can disagree with the cut without having to reverse-engineer it.
#
# What this is NOT: a measurement of contamination. decontX estimates that with an
# EM model over a chosen ambient profile. This is a cheap screen that says "worth
# looking at", and the action's evidence node is what carries the estimate.
sclet_ai_ambient_pressure <- function(object) {
    threshold_umi_q10 <- 0.05
    threshold_detection <- 0.90

    counts <- tryCatch(SummarizedExperiment::assay(object, "counts"),
        error = function(e) NULL)
    if (is.null(counts)) {
        return(list(
            available = FALSE,
            reason = "counts_assay_missing",
            note = "ambient pressure is screened from raw counts",
            raw_values_included = FALSE
        ))
    }

    # Per-cell total counts, as a plain vector. Sparse summary() is used for the
    # column sums rather than converting the matrix.
    totals <- tryCatch({
        if (inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix")) {
            Matrix::colSums(counts)
        } else {
            colSums(as.matrix(counts))
        }
    }, error = function(e) NULL)
    if (is.null(totals) || !length(totals)) {
        return(list(available = FALSE, reason = "counts_unreadable",
            raw_values_included = FALSE))
    }

    # D1: how much of the library sits in the low-UMI tail. Ambient RNA is
    # proportionally louder in cells with few real transcripts, so a heavy lower
    # tail is where contamination shows first.
    q10 <- stats::quantile(totals, 0.10, na.rm = TRUE, names = FALSE)
    low_tail_share <- if (sum(totals) > 0) {
        sum(totals[totals <= q10]) / sum(totals)
    } else {
        NA_real_
    }

    # D2: among the most abundant genes, how uniformly are they detected? A gene
    # that appears in almost every cell regardless of cell type is the signature
    # of a shared ambient pool (or a housekeeping gene -- this screen cannot tell
    # those apart, and says so).
    gene_totals <- tryCatch({
        if (inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix")) {
            Matrix::rowSums(counts)
        } else {
            rowSums(as.matrix(counts))
        }
    }, error = function(e) NULL)
    detection_uniformity <- NA_real_
    if (!is.null(gene_totals) && length(gene_totals)) {
        top <- order(gene_totals, decreasing = TRUE)
        top <- utils::head(top, 200L)
        if (length(top)) {
            sub <- counts[top, , drop = FALSE]
            detected <- tryCatch({
                if (inherits(sub, "dgCMatrix") || inherits(sub, "dgTMatrix")) {
                    Matrix::rowSums(sub > 0) / ncol(sub)
                } else {
                    rowMeans(as.matrix(sub) > 0)
                }
            }, error = function(e) NULL)
            if (!is.null(detected) && length(detected)) {
                detection_uniformity <- mean(detected, na.rm = TRUE)
            }
        }
    }

    # A coarse label from explicit cuts. `high` requires BOTH screens to fire:
    # a heavy low-UMI tail alone is normal for many tissues, and a high detection
    # rate alone is what housekeeping genes look like.
    low_tail_fires <- !is.na(low_tail_share) && low_tail_share <= threshold_umi_q10
    uniform_fires <- !is.na(detection_uniformity) &&
        detection_uniformity >= threshold_detection
    label <- if (isTRUE(low_tail_fires) && isTRUE(uniform_fires)) {
        "high"
    } else if (isTRUE(uniform_fires) || isTRUE(low_tail_fires)) {
        "moderate"
    } else {
        "low"
    }

    list(
        available = TRUE,
        # the two screens, each with its own name and its own meaning
        low_umi_tail_share = if (is.na(low_tail_share)) NA_real_ else round(low_tail_share, 4L),
        top_gene_detection_uniformity = if (is.na(detection_uniformity)) {
            NA_real_
        } else {
            round(detection_uniformity, 4L)
        },
        # and the cuts that turn them into a label
        threshold_umi_q10 = threshold_umi_q10,
        threshold_detection = threshold_detection,
        label = label,
        caveat = paste(
            "This is a screen, not a measurement: a high detection rate is also",
            "what housekeeping genes look like, and a heavy low-UMI tail is",
            "normal for many tissues. 'high' means both screens fired and",
            "cleaning is worth considering -- the contamination estimate itself",
            "comes from the run_decontx evidence node, not from here."
        ),
        raw_values_included = FALSE
    )
}
