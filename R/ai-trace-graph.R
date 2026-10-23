# The execution trace as a graph.
#
# WHY THIS FILE EXISTS: the plan is a DAG (`depends_on`, `${step_6.output...}`
# references), but until this change the execution record was a flat array and
# nothing linked a step to the records it wrote. A finished round could be read
# back as a list of steps and a list of results, never as a path. This is a
# read-time projection over the records that already exist -- no second copy of
# anything, no storage rewrite -- plus the queries that make it usable.
#
# WHAT IT CANNOT DO, on purpose: both record stores overwrite by id
# (R/state.R:824, :633), so this is a current-state DAG plus the trace, not a
# version history. A re-run with the same id shows two producers for one node;
# the node holds the latest values.

sclet_ai_graph_chr <- function(value) {
    if (is.null(value) || !length(value)) {
        return(NA_character_)
    }
    out <- as.character(value)[[1L]]
    if (is.na(out)) NA_character_ else out
}

sclet_ai_graph_time <- function(value) {
    sclet_ai_graph_chr(value)
}

# `sclet_ai_graph_chr()` returns NA for a missing value, which `%||%` does not
# replace (NA is not NULL). Labels and types need the fallback form.
sclet_ai_graph_chr_or <- function(value, fallback) {
    out <- sclet_ai_graph_chr(value)
    if (is.na(out) || !nzchar(out)) fallback else out
}

sclet_ai_graph_number <- function(value) {
    if (is.null(value) || !length(value)) {
        return(NA_real_)
    }
    out <- suppressWarnings(as.numeric(value)[[1L]])
    if (is.na(out)) NA_real_ else out
}

sclet_ai_graph_integer <- function(value) {
    if (is.null(value) || !length(value)) {
        return(NA_integer_)
    }
    out <- suppressWarnings(as.integer(value)[[1L]])
    if (is.na(out)) NA_integer_ else out
}

sclet_ai_graph_node_id <- function(type, id) {
    if (identical(as.character(type), "ai_evidence")) {
        return(paste0("evidence:", id))
    }
    paste0("record:", type, "::", id)
}

# Map a key as stored in a record (`"<type>::<id>"`, `"analysis::<id>"`) to a
# node id. Evidence is always the `evidence:` node, never also a `record:` node:
# one mapping function, so a step edge and an evidence edge cannot disagree.
sclet_ai_graph_key_node <- function(key) {
    key <- as.character(key)
    parts <- strsplit(key, "::", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) {
        return(NA_character_)
    }
    type <- parts[[1L]]
    id <- paste(parts[-1L], collapse = "::")
    if (identical(type, "ai_evidence")) {
        return(paste0("evidence:", id))
    }
    if (identical(type, "analysis")) {
        return(paste0("record:analysis::", id))
    }
    paste0("record:", type, "::", id)
}

sclet_ai_graph_node_row <- function(id, kind, label, type, status, created_at,
                                    complete = NA, action = NA_character_,
                                    duration_sec = NA_real_,
                                    attempts = NA_integer_) {
    list(
        id = as.character(id),
        kind = as.character(kind),
        label = as.character(label),
        type = as.character(type),
        status = as.character(status),
        created_at = as.character(created_at),
        complete = complete,
        action = as.character(action),
        duration_sec = as.numeric(duration_sec),
        attempts = as.integer(attempts)
    )
}

sclet_ai_graph_edge_row <- function(from, to, kind,
                                    parent_kind = NA_character_,
                                    declared_input = NA_character_) {
    list(
        from = as.character(from),
        to = as.character(to),
        kind = as.character(kind),
        parent_kind = as.character(parent_kind),
        declared_input = as.character(declared_input)
    )
}

# Flatten character leaves of a record's `inputs`, keeping the field path.
sclet_ai_graph_input_values <- function(value, path = character(), depth = 0L) {
    if (depth > 3L || is.null(value) || !length(value)) {
        return(list())
    }
    if (is.list(value)) {
        names_value <- names(value)
        out <- list()
        for (i in seq_along(value)) {
            field <- if (!is.null(names_value) && nzchar(names_value[[i]])) {
                names_value[[i]]
            } else {
                as.character(i)
            }
            out <- c(out, sclet_ai_graph_input_values(value[[i]], c(path, field), depth + 1L))
        }
        return(out)
    }
    if (is.character(value) && length(value) <= 8L) {
        return(list(list(path = paste(path, collapse = "."), values = value)))
    }
    list()
}

# Fields that hold provenance bookkeeping, not record references. A report's
# `inputs$task` is "report"; a record named "report" would otherwise be cited as
# an input by every report on the object.
sclet_ai_graph_input_is_metadata <- function(path) {
    path <- tolower(as.character(path))
    last <- sub(".*\\.", "", path)
    last %in% c("task", "plan_id", "action_ids", "context_fingerprint",
                "context_schema_version") ||
        grepl("fingerprint|schema_version", path)
}

# Prefer the type a field name declares: `detest_id = "markers"` is a detest
# record, even when a graph record happens to be called "markers" too.
sclet_ai_graph_typed_target <- function(nodes, path, value) {
    last <- tolower(sub(".*\\.", "", as.character(path)))
    if (!grepl("_id$", last)) {
        return(NA_character_)
    }
    type <- sub("_id$", "", last)
    if (!nzchar(type)) {
        return(NA_character_)
    }
    candidate <- paste0("record:", type, "::", as.character(value))
    if (nzchar(candidate) && !is.null(nodes[[candidate]])) candidate else NA_character_
}

sclet_ai_graph_frame <- function(rows, template) {
    if (!length(rows)) {
        frame <- as.data.frame(template, stringsAsFactors = FALSE)
        return(frame[0L, , drop = FALSE])
    }
    frames <- lapply(rows, function(row) as.data.frame(row, stringsAsFactors = FALSE))
    frame <- do.call(rbind, frames)
    rownames(frame) <- NULL
    frame[, names(template), drop = FALSE]
}

sclet_ai_graph_build <- function(object) {
    state <- sclet_get_state(object)
    records <- state$states$records
    if (!is.list(records)) {
        records <- list()
    }
    analyses <- state$analyses
    if (!is.list(analyses)) {
        analyses <- list()
    }

    nodes <- list()
    edges <- list()
    node_template <- sclet_ai_graph_node_row(
        NA_character_, NA_character_, NA_character_, NA_character_,
        NA_character_, NA_character_
    )
    edge_template <- sclet_ai_graph_edge_row(
        NA_character_, NA_character_, NA_character_
    )
    add_node <- function(row) {
        nodes[[row$id]] <<- row
        invisible(NULL)
    }
    add_edge <- function(row) {
        edges[[length(edges) + 1L]] <<- row
        invisible(NULL)
    }

    bare_index <- list()
    analysis_index <- list()
    index_bare <- function(id, node_id) {
        id <- as.character(id)
        if (!nzchar(id)) {
            return(invisible(NULL))
        }
        current <- bare_index[[id]]
        if (is.null(current)) {
            bare_index[[id]] <<- node_id
        } else if (!identical(current, node_id)) {
            bare_index[[id]] <<- "" # ambiguous: two record kinds share this id
        }
        invisible(NULL)
    }

    # ---- state records: analysis results and evidence --------------------
    for (type in names(records)) {
        type_records <- records[[type]]
        if (!is.list(type_records) || !length(type_records)) {
            next
        }
        ids <- names(type_records)
        if (is.null(ids)) {
            next
        }
        for (id in ids[nzchar(ids)]) {
            record <- type_records[[id]]
            is_evidence <- identical(type, "ai_evidence")
            node_id <- sclet_ai_graph_node_id(type, id)
            add_node(sclet_ai_graph_node_row(
                id = node_id,
                kind = if (is_evidence) "evidence" else "record",
                label = if (is_evidence) id else paste0(type, "::", id),
                type = type,
                status = if (is.list(record)) sclet_ai_graph_chr(record$status) else NA_character_,
                created_at = if (is.list(record)) sclet_ai_graph_time(record$created_at) else NA_character_
            ))
            index_bare(id, node_id)
        }
    }

    # ---- analysis store: execution rounds, and everything else -----------
    execution_ids <- character()
    other_analysis_ids <- character()
    for (id in names(analyses)) {
        record <- analyses[[id]]
        if (!is.list(record)) {
            next
        }
        record_type <- sclet_ai_graph_chr(record$type)
        if (identical(record_type, "ai_execution")) {
            execution_ids <- c(execution_ids, id)
        } else {
            other_analysis_ids <- c(other_analysis_ids, id)
        }
    }
    for (id in other_analysis_ids) {
        record <- analyses[[id]]
        node_id <- paste0("record:analysis::", id)
        add_node(sclet_ai_graph_node_row(
            id = node_id,
            kind = "record",
            label = paste0("analysis::", id),
            type = sclet_ai_graph_chr_or(record$type, "analysis"),
            status = sclet_ai_graph_chr(record$summary$status %||% record$status),
            created_at = sclet_ai_graph_time(record$created_at)
        ))
        index_bare(id, node_id)
        analysis_index[[id]] <- node_id
    }

    # ---- rounds and steps ------------------------------------------------
    unresolved <- 0L
    legacy_rounds <- 0L
    object_replaced <- FALSE
    round_stamps <- vapply(execution_ids, function(id) {
        sclet_ai_graph_chr_or(analyses[[id]]$created_at, "")
    }, character(1L))
    execution_ids <- execution_ids[order(round_stamps, execution_ids)]
    for (id in execution_ids) {
        record <- analyses[[id]]
        round_id <- paste0("round:", id)
        steps <- record$params$steps %||% list()
        # A round with no recorded steps has no topology to trust: an execution
        # record can exist after a plan that never ran a step. Claiming TRUE
        # here would say the trace was recorded when there is nothing to record.
        # The same holds for a store that was written by something other than
        # this package version: a malformed `params$steps` must degrade to
        # "incomplete", never to an error in a read-only projection.
        elements_ok <- is.list(steps) &&
            all(vapply(steps, is.list, logical(1L)))
        if (!is.list(steps)) {
            steps <- list()
        }
        # A store written by something else can hold a non-list element next to
        # real steps; keep the real ones, and mark the round incomplete rather
        # than erroring in a read-only projection.
        steps <- Filter(function(step) is.list(step), steps)
        steps_ok <- length(steps) > 0L && elements_ok &&
            all(vapply(steps, function(step) {
                nzchar(as.character(step$id %||% ""))
            }, logical(1L)))
        complete <- steps_ok &&
            all(vapply(steps, function(step) !is.null(step$depends_on), logical(1L)))
        if (!complete) {
            legacy_rounds <- legacy_rounds + 1L
        }
        add_node(sclet_ai_graph_node_row(
            id = round_id,
            kind = "round",
            label = sclet_ai_graph_chr_or(record$inputs$plan_id, id),
            type = "ai_execution",
            status = sclet_ai_graph_chr(record$summary$status %||% record$status),
            created_at = sclet_ai_graph_time(record$created_at),
            complete = complete
        ))
        index_bare(id, round_id)
        analysis_index[[id]] <- round_id
        step_ids <- vapply(steps, function(step) as.character(step$id %||% ""), character(1L))
        for (step in steps) {
            step_id <- as.character(step$id %||% "")
            if (!nzchar(step_id)) {
                next
            }
            node_id <- paste0("step:", id, ":", step_id)
            add_node(sclet_ai_graph_node_row(
                id = node_id,
                kind = "step",
                label = sclet_ai_graph_chr_or(step$action, step_id),
                type = "step",
                status = sclet_ai_graph_chr(step$status),
                created_at = sclet_ai_graph_time(record$created_at),
                complete = complete,
                action = sclet_ai_graph_chr(step$action),
                duration_sec = sclet_ai_graph_number(step$duration_sec),
                attempts = sclet_ai_graph_integer(step$attempts)
            ))
            if (isTRUE(step$replaces_object)) {
                object_replaced <- TRUE
            }
            add_edge(sclet_ai_graph_edge_row(node_id, round_id, "step_of"))
            for (dependency in as.character(step$depends_on %||% character())) {
                target <- paste0("step:", id, ":", dependency)
                if (nzchar(dependency) && dependency %in% step_ids) {
                    add_edge(sclet_ai_graph_edge_row(node_id, target, "depends_on"))
                } else {
                    unresolved <- unresolved + 1L
                }
            }
            for (key in as.character(step$record_ids %||% character())) {
                target <- sclet_ai_graph_key_node(key)
                if (!is.na(target) && !is.null(nodes[[target]])) {
                    add_edge(sclet_ai_graph_edge_row(node_id, target, "produced"))
                } else {
                    unresolved <- unresolved + 1L
                }
            }
            for (key in as.character(step$refreshed_ids %||% character())) {
                target <- sclet_ai_graph_key_node(key)
                if (!is.na(target) && !is.null(nodes[[target]])) {
                    add_edge(sclet_ai_graph_edge_row(node_id, target, "refreshed"))
                } else {
                    unresolved <- unresolved + 1L
                }
            }
        }
    }

    resolve_bare <- function(id) {
        if (is.null(id) || !length(id)) {
            return(NA_character_)
        }
        id <- as.character(id)[[1L]]
        if (is.na(id) || !nzchar(id)) {
            return(NA_character_)
        }
        target <- bare_index[[id]]
        if (is.null(target) || !nzchar(target)) {
            return(NA_character_)
        }
        target
    }

    # An evidence node's `source` names the analysis that produced it. A state
    # record and an analysis record can share a bare id (RunSuperCell writes
    # "supercell" in both stores), which makes the bare index ambiguous; for a
    # source, the analysis store is the answer.
    resolve_source <- function(id) {
        if (is.null(id) || !length(id)) {
            return(NA_character_)
        }
        id <- as.character(id)[[1L]]
        if (is.na(id) || !nzchar(id)) {
            return(NA_character_)
        }
        target <- analysis_index[[id]]
        if (!is.null(target) && nzchar(target)) {
            return(target)
        }
        resolve_bare(id)
    }

    # ---- record -> record: declared inputs -------------------------------
    for (type in names(records)) {
        # An evidence node's `inputs` mirror `source` and `parents`, which
        # `source_of` / `derived_from` already express with a precise kind. A
        # `consumed` edge here would duplicate them under a vaguer name.
        if (identical(type, "ai_evidence")) {
            next
        }
        type_records <- records[[type]]
        if (!is.list(type_records)) {
            next
        }
        for (id in names(type_records)) {
            record <- type_records[[id]]
            if (!is.list(record) || !is.list(record$inputs)) {
                next
            }
            from <- sclet_ai_graph_node_id(type, id)
            for (leaf in sclet_ai_graph_input_values(record$inputs)) {
                if (sclet_ai_graph_input_is_metadata(leaf$path)) {
                    next
                }
                for (value in leaf$values) {
                    target <- sclet_ai_graph_typed_target(nodes, leaf$path, value)
                    if (is.na(target)) {
                        target <- resolve_bare(value)
                    }
                    if (is.na(target) || identical(target, from)) {
                        next
                    }
                    add_edge(sclet_ai_graph_edge_row(
                        from, target, "consumed",
                        declared_input = leaf$path
                    ))
                }
            }
        }
    }
    for (id in other_analysis_ids) {
        record <- analyses[[id]]
        if (!is.list(record$inputs)) {
            next
        }
        from <- paste0("record:analysis::", id)
        for (leaf in sclet_ai_graph_input_values(record$inputs)) {
            if (sclet_ai_graph_input_is_metadata(leaf$path)) {
                next
            }
            for (value in leaf$values) {
                target <- sclet_ai_graph_typed_target(nodes, leaf$path, value)
                if (is.na(target)) {
                    target <- resolve_bare(value)
                }
                if (is.na(target) || identical(target, from)) {
                    next
                }
                add_edge(sclet_ai_graph_edge_row(
                    from, target, "consumed",
                    declared_input = leaf$path
                ))
            }
        }
    }

    # ---- evidence -> source / parents ------------------------------------
    evidence <- records$ai_evidence
    if (is.list(evidence) && length(evidence)) {
        for (id in names(evidence)) {
            record <- evidence[[id]]
            if (!is.list(record)) {
                next
            }
            # An evidence node stores its own fields in `summary` (written by
            # RecordAIEvidence through sclet_set_analysis_state), with `source`
            # and `parents` also mirrored in `inputs`. Read the summary first so
            # a hand-built or legacy node that keeps them at the top level still
            # resolves.
            node <- if (is.list(record$summary)) record$summary else record
            from <- paste0("evidence:", id)
            # A source that resolves is drawn, every one of them; an empty or
            # missing source names nothing and is not a defect.
            sources <- node$source %||% record$inputs$source
            sources <- if (is.null(sources)) character() else
                as.character(unlist(sources, use.names = FALSE))
            sources <- unique(sources[!is.na(sources) & nzchar(sources)])
            seen_sources <- character()
            for (source_value in sources) {
                source <- resolve_source(source_value)
                if (is.na(source)) {
                    unresolved <- unresolved + 1L
                    next
                }
                if (source %in% seen_sources) {
                    next
                }
                seen_sources <- c(seen_sources, source)
                add_edge(sclet_ai_graph_edge_row(
                    from, source, "source_of",
                    parent_kind = nodes[[source]]$kind
                ))
            }
            # A parent that is listed but cannot resolve is a real defect, so
            # it counts; a repeated parent counts once.
            raw_parents <- node$parents %||% record$inputs$parents %||% character()
            if (!is.atomic(raw_parents)) {
                raw_parents <- unlist(raw_parents, use.names = FALSE)
            }
            if (is.null(raw_parents)) {
                raw_parents <- character()
            }
            parents <- unique(as.character(raw_parents))
            seen_parents <- character()
            for (parent in parents) {
                target <- resolve_bare(parent)
                if (is.na(target)) {
                    unresolved <- unresolved + 1L
                    next
                }
                if (target %in% seen_parents) {
                    next
                }
                seen_parents <- c(seen_parents, target)
                add_edge(sclet_ai_graph_edge_row(
                    from, target, "derived_from",
                    parent_kind = nodes[[target]]$kind
                ))
            }
        }
    }

    structure(
        list(
            meta = list(
                schema = "sclet-trace-graph/1",
                n_rounds = length(execution_ids),
                legacy_rounds = legacy_rounds,
                unresolved_edges = unresolved,
                # TRUE when a recorded step returned a different object (a
                # metacell object, for example). Earlier rounds are NOT lost by
                # this projection -- they are gone from the object, and this flag
                # is how a reader learns that the absence is not a projection bug.
                object_replaced = object_replaced
            ),
            nodes = sclet_ai_graph_frame(nodes, node_template),
            edges = sclet_ai_graph_frame(edges, edge_template)
        ),
        class = c("sclet_trace_graph", "list")
    )
}

#' The execution trace as a graph
#'
#' Project the recorded AI rounds, their steps, the analysis records they wrote
#' and the evidence nodes onto one typed node/edge table.
#'
#' The projection is read-only and is built from the records the object already
#' carries: no second copy, no storage rewrite. `depends_on` edges come from
#' what `ExecuteAIPlan()` stored for each step, `produced` / `refreshed` edges
#' from the record diff taken at the step's commit point, `consumed` edges only
#' from ids a record already declares in its `inputs`, and `source_of` /
#' `derived_from` from the evidence nodes themselves. `source_of` targets a
#' `round` node when an evidence node's source names an execution record rather
#' than an analysis record; `parent_kind` says which.
#'
#' Both record stores overwrite by id, so this is a current-state DAG plus the
#' trace, not a version history: a re-run with the same id shows two producers
#' for one node, and the node holds the latest values. A round recorded before
#' the trace fields existed has `complete = FALSE`; its evidence edges are still
#' present and the queries report the incompleteness.
#'
#' The tables are a structural whitelist, not a sanitised copy: the ten columns
#' are drawn from record ids, record types, statuses, action names and the plan
#' id, and nothing is read from `artifacts` or from any analysis value, matrix
#' or cell-level string. `status` and the round label are copied verbatim from
#' record metadata, so a caller that writes free text into a record's `status`
#' would see it here; if you need a recorded value, read the record.
#'
#' @param object A `SingleCellExperiment` with at least one recorded AI round.
#' @return A `sclet_trace_graph`: a list with `meta` (schema, round counts and
#'   `unresolved_edges`, the number of references that could not be drawn
#'   because their target did not exist: a dangling step dependency, record id,
#'   evidence source or evidence parent, and `object_replaced`, TRUE when a
#'   recorded step returned a different object and the earlier records are
#'   therefore gone rather than unprojected), `nodes` (one row per round, step,
#'   record and evidence
#'   node) and `edges` (one row per edge, with `kind`).
#' @seealso [graph_neighbors()], [graph_producers()], [graph_consumers()]
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     list(counts = SummarizedExperiment::assay(sce, "counts"))
#' )
#' registry <- AIDefaultExecutionRegistry(sce, include = "read")
#' plan <- list(
#'     plan_id = "plan_trace_demo",
#'     task = "inspect",
#'     actions = list(list(id = "s1", action = "inspect_status", params = list()))
#' )
#' validation <- ValidateAIPlan(plan, object = sce, registry = registry)
#' executed <- ExecuteAIPlan(
#'     sce, plan, registry,
#'     validation = validation, dry_run = FALSE,
#'     confirmation = validation$confirmation_token
#' )
#' graph <- AnalysisGraph(executed$object)
#' graph$meta
#' graph$edges
#' @export
AnalysisGraph <- function(object) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.", call. = FALSE)
    }
    sclet_ai_graph_build(object)
}

sclet_ai_graph_check <- function(graph) {
    if (!inherits(graph, "sclet_trace_graph")) {
        stop("`graph` must come from AnalysisGraph().", call. = FALSE)
    }
    invisible(graph)
}

# Are the rounds behind these nodes fully recorded?
sclet_ai_graph_complete <- function(graph, ids) {
    nodes <- graph$nodes[graph$nodes$id %in% ids, , drop = FALSE]
    rounds <- nodes$id[nodes$kind == "round"]
    steps <- nodes$id[nodes$kind == "step"]
    if (length(steps)) {
        # A step id is `step:<execution>:<step>`; its round node id is
        # `round:<execution>`, and the prefix is what makes the lookup match.
        rounds <- c(rounds, paste0("round:", sub("^step:([^:]+):.*$", "\\1", steps)))
    }
    rounds <- unique(rounds[nzchar(rounds)])
    if (!length(rounds)) {
        # A record with no producer edge involves no round. The honest answer
        # is how completely the graph as a whole is recorded, not "TRUE": on a
        # legacy object that answer would claim a trace that does not exist.
        rounds <- graph$nodes$id[graph$nodes$kind == "round"]
    }
    if (!length(rounds)) {
        # No round at all is "no recorded trace", not "nothing missing".
        return(FALSE)
    }
    flags <- graph$nodes$complete[match(rounds, graph$nodes$id)]
    all(!is.na(flags) & flags)
}

#' Neighbours of one trace-graph node
#'
#' @param graph A graph returned by `AnalysisGraph()`.
#' @param id A node id, e.g. `"step:<execution>:s2"`.
#' @param direction `"out"` (edges from `id`), `"in"` (edges into `id`) or
#'   `"both"`.
#' @return A list with `ids` (neighbour node ids) and `complete` (whether every
#'   round involved carries the full trace).
#' @examples
#' # see AnalysisGraph()
#' @export
graph_neighbors <- function(graph, id, direction = c("out", "in", "both")) {
    sclet_ai_graph_check(graph)
    direction <- match.arg(direction)
    edges <- graph$edges
    id <- as.character(id)
    out <- character()
    if (direction %in% c("out", "both")) {
        out <- c(out, edges$to[edges$from %in% id])
    }
    if (direction %in% c("in", "both")) {
        out <- c(out, edges$from[edges$to %in% id])
    }
    out <- unique(out[nzchar(out)])
    list(ids = out, complete = sclet_ai_graph_complete(graph, c(id, out)))
}

#' Which steps produced a record
#'
#' @param graph A graph returned by `AnalysisGraph()`.
#' @param record_id A record or evidence node id, e.g. `"record:detest::markers"`.
#' @return A list with `ids` (step node ids) and `complete`.
#' @examples
#' # see AnalysisGraph()
#' @export
graph_producers <- function(graph, record_id) {
    sclet_ai_graph_check(graph)
    edges <- graph$edges
    produced <- edges$kind %in% c("produced", "refreshed")
    ids <- unique(edges$from[produced & edges$to %in% as.character(record_id)])
    list(ids = ids, complete = sclet_ai_graph_complete(graph, c(as.character(record_id), ids)))
}

#' Which later steps declared a step as upstream
#'
#' @param graph A graph returned by `AnalysisGraph()`.
#' @param step_id A step node id, e.g. `"step:<execution>:s2"`.
#' @return A list with `ids` (step node ids that depend on `step_id`) and
#'   `complete`.
#' @examples
#' # see AnalysisGraph()
#' @export
graph_consumers <- function(graph, step_id) {
    sclet_ai_graph_check(graph)
    edges <- graph$edges
    ids <- unique(edges$from[edges$kind %in% "depends_on" & edges$to %in% as.character(step_id)])
    list(ids = ids, complete = sclet_ai_graph_complete(graph, c(as.character(step_id), ids)))
}

#' @export
print.sclet_trace_graph <- function(x, ...) {
    meta <- x$meta
    cat("<sclet trace graph>", meta$schema, "\n")
    cat("  rounds:", meta$n_rounds,
        if (meta$legacy_rounds > 0L) paste0(" (", meta$legacy_rounds, " without stored topology)") else "",
        "\n")
    cat("  nodes:", nrow(x$nodes), "| edges:", nrow(x$edges), "\n")
    if (meta$unresolved_edges > 0L) {
        cat("  unresolved edges (dangling target):", meta$unresolved_edges, "\n")
    }
    if (isTRUE(meta$object_replaced)) {
        cat("  note: a step replaced the object; earlier rounds are not on it\n")
    }
    invisible(x)
}
