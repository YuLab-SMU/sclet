sclet_ai_claim_levels <- c("observed", "measured", "associated", "consistent_with")

sclet_ai_claim_rank <- function(level) {
    match(as.character(level), sclet_ai_claim_levels)
}

sclet_ai_claim_from_rank <- function(rank) {
    if (is.na(rank) || rank < 1L) return(NULL)
    sclet_ai_claim_levels[[min(as.integer(rank), length(sclet_ai_claim_levels))]]
}

sclet_ai_claim_ancestors <- function(all_nodes, id) {
    seen <- character()
    queue <- id
    while (length(queue)) {
        current <- queue[[1L]]
        queue <- queue[-1L]
        if (current %in% seen) next
        seen <- c(seen, current)
        node <- all_nodes[[current]]
        parents <- as.character(node$parents %||% character())
        parents <- parents[parents %in% names(all_nodes)]
        queue <- c(queue, setdiff(parents, seen))
    }
    seen
}

sclet_ai_claim_dependency_reasons <- function(all_nodes, a, b) {
    node_a <- all_nodes[[a]]
    node_b <- all_nodes[[b]]
    reasons <- character()
    group_a <- as.character(node_a$dependency_group %||% "")
    group_b <- as.character(node_b$dependency_group %||% "")
    if (nzchar(group_a) && identical(group_a, group_b)) {
        reasons <- c(reasons, "shared_dependency_group")
    }
    source_a <- as.character(node_a$source %||% "")
    source_b <- as.character(node_b$source %||% "")
    # two nodes produced by the same analysis execution are one measurement,
    # even when their dependency_group hashes differ
    if (nzchar(source_a) && identical(source_a, source_b)) {
        reasons <- c(reasons, "shared_source")
    }
    if (b %in% setdiff(sclet_ai_claim_ancestors(all_nodes, a), a)) {
        reasons <- c(reasons, "ancestor_chain")
    }
    if (a %in% setdiff(sclet_ai_claim_ancestors(all_nodes, b), b)) {
        reasons <- c(reasons, "descendant_chain")
    }
    reasons
}
#' Audit whether a claim is supported by independent evidence
#'
#' `sclet_ai_claim_ceiling()` answers a narrow question: given the evidence nodes
#' recorded on an object, what is the strongest claim those nodes can actually
#' carry? It never produces or strengthens a claim. It groups the supplied nodes
#' into lines of independent support and returns the ceiling implied by that
#' structure.
#'
#' Two nodes are treated as dependent - and therefore as one line of support -
#' when they share a `dependency_group`, when one descends from the other through
#' the `parents` chain, or when both were produced by the same underlying
#' analysis (the same `source`). The last rule matters: a single run can emit
#' many evidence nodes, and counting them separately would inflate apparent
#' support. Groups are then formed as connected components of that relation,
#' which is the conservative direction - it can never overstate independence.
#'
#' The ceiling is derived deterministically. With no line of support nothing
#' may be claimed. A single line of support is capped at `"associated"` so a
#' lone node can never be upgraded into a stronger statement. With two or more
#' independent lines the ceiling is at least `"associated"` (a synthesis of
#' independent sources is itself an inference, never a direct measurement) and
#' at most `"consistent_with"`, and it is additionally bounded by the weakest
#' node in the set. `user_decision` nodes are excluded from the support count,
#' matching the existing rule that an AI-generated claim may not cite a human
#' decision as evidence.
#'
#' @param object A `SingleCellExperiment` object.
#' @param evidence_ids Character vector of evidence ids to audit. Defaults to
#'   every completed evidence node on the object.
#' @param proposed_claim_level Optional claim level to check against the
#'   ceiling, one of `"observed"`, `"measured"`, `"associated"` or
#'   `"consistent_with"`.
#' @return A list with the `ceiling_claim_level`, the number of independent
#'   `n_independent_groups`, the supporting `groups`, a `dependent_pairs`
#'   diagnostic, and a `status` of `"unsupported"`, `"allowed"` or
#'   `"downgraded"`. This function adds no evidence and modifies nothing.
#' @noRd
sclet_ai_claim_ceiling <- function(object, evidence_ids = NULL, proposed_claim_level = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    if (!is.null(proposed_claim_level) &&
        (!is.character(proposed_claim_level) || length(proposed_claim_level) != 1L ||
            !proposed_claim_level %in% sclet_ai_claim_levels)) {
        stop("proposed_claim_level must be one of: ", paste(sclet_ai_claim_levels, collapse = ", "))
    }
    empty <- function(reason, n_nodes = 0L) {
        list(
            status = "unsupported",
            ceiling_claim_level = NULL,
            requested_claim_level = proposed_claim_level,
            n_nodes = n_nodes,
            n_independent_groups = 0L,
            groups = list(),
            dependent_pairs = list(),
            reasons = reason,
            raw_values_included = FALSE
        )
    }
    all_nodes <- sclet_ai_evidence_get_all(object)
    if (!length(all_nodes)) return(empty("no evidence has been recorded on this object"))
    completed <- vapply(all_nodes, function(node) {
        identical(as.character(node$status %||% "completed"), "completed")
    }, logical(1L))
    ids <- names(all_nodes)[completed]
    if (!is.null(evidence_ids)) {
        missing <- setdiff(as.character(evidence_ids), ids)
        if (length(missing)) {
            stop("unknown or incomplete evidence id(s): ", paste(missing, collapse = ", "))
        }
        ids <- intersect(as.character(evidence_ids), ids)
    }
    if (!length(ids)) return(empty("no completed evidence node was supplied"))
    nodes <- all_nodes[ids]

    user_decided <- ids[vapply(nodes, function(node) {
        sclet_ai_evidence_kind_is_user_decision(node$kind %||% "")
    }, logical(1L))]
    supporting <- setdiff(ids, user_decided)

    dependent_pairs <- list()
    if (length(ids) > 1L) {
        for (i in seq_len(length(ids) - 1L)) {
            for (j in seq(i + 1L, length(ids))) {
                reasons <- sclet_ai_claim_dependency_reasons(all_nodes, ids[[i]], ids[[j]])
                if (length(reasons)) {
                    dependent_pairs[[length(dependent_pairs) + 1L]] <- list(
                        a = ids[[i]], b = ids[[j]], reasons = reasons
                    )
                }
            }
        }
    }

# connected components over the dependency relation
    group_of <- stats::setNames(rep(NA_integer_, length(supporting)), supporting)
    n_groups <- 0L
    for (id in supporting) {
        if (!is.na(group_of[[id]])) next
        n_groups <- n_groups + 1L
        group_of[[id]] <- n_groups
        queue <- id
        while (length(queue)) {
            current <- queue[[1L]]
            queue <- queue[-1L]
            for (other in supporting) {
                if (!is.na(group_of[[other]])) next
                if (length(sclet_ai_claim_dependency_reasons(all_nodes, current, other))) {
                    group_of[[other]] <- n_groups
                    queue <- c(queue, other)
                }
            }
        }
    }
    groups <- lapply(seq_len(n_groups), function(index) {
        members <- names(group_of)[group_of == index]
        list(
            group = index,
            n_nodes = length(members),
            nodes = members,
            sources = sort(unique(vapply(nodes[members], function(n) as.character(n$source %||% ""), character(1L))))
        )
    })

    if (!length(supporting)) {
        result <- empty(
            "only user_decision evidence was supplied; a human decision cannot support an AI-generated claim",
            n_nodes = length(ids)
        )
        result$dependent_pairs <- dependent_pairs
        result$user_decision_nodes <- user_decided
        return(result)
    }

    ranks <- vapply(nodes[supporting], function(n) {
        sclet_ai_claim_rank(n$claim_level %||% "observed")
    }, integer(1L))
    min_rank <- min(ranks)
    if (n_groups <= 1L) {
        # a single line of support may never be upgraded into a stronger claim
        ceiling_rank <- min(3L, min_rank)
    } else {
        # a synthesis of independent sources is an inference, never a measurement
        ceiling_rank <- max(3L, min_rank)
    }
    ceiling <- sclet_ai_claim_from_rank(ceiling_rank)

    reasons <- character()
    if (n_groups <= 1L) {
        reasons <- c(reasons, paste0(
            "all supporting evidence resolves to ", n_groups,
            " independent line(s) of support; a claim cannot exceed 'associated'"
        ))
    } else {
        reasons <- c(reasons, paste0(
            n_groups, " independent lines of support; the ceiling is additionally",
            " bounded by the weakest supporting node"
        ))
    }
    if (length(user_decided)) {
        reasons <- c(reasons, "user_decision nodes were excluded from the support count")
    }

    status <- "allowed"
    if (!is.null(proposed_claim_level)) {
        if (sclet_ai_claim_rank(proposed_claim_level) > ceiling_rank) {
            status <- "downgraded"
            reasons <- c(reasons, paste0(
                "requested '", proposed_claim_level, "' exceeds the ceiling '", ceiling,
                "'; use '", ceiling, "' or weaker"
            ))
        }
    }
    list(
        status = status,
        ceiling_claim_level = ceiling,
        requested_claim_level = proposed_claim_level,
        n_nodes = length(ids),
        n_independent_groups = n_groups,
        groups = groups,
        dependent_pairs = dependent_pairs,
        user_decision_nodes = user_decided,
        reasons = reasons,
        raw_values_included = FALSE
    )
}