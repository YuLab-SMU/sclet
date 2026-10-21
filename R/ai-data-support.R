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
