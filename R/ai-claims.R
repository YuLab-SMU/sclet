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
#' Finding claim levels that must cite independent support
#'
#' A finding at `"hypothesis"` is explicitly speculative and self-limiting, so
#' it needs no evidence to be recorded: it asserts nothing. The levels that
#' actually assert something beyond one node's own claim level are the ones
#' gated here, plus `"causal"`, which is banned outright below.
sclet_ai_findings_needing_support <- c("consistent_with", "suggestive", "causal")

# Map a finding's claim level onto the strongest evidence ceiling it may cite.
# "suggestive" sits one step above "consistent_with" on the finding ladder, so
# it must clear the same bar. "causal" never reaches here: it is rejected
# outright above.
sclet_ai_finding_allowed_level <- function(claim_level) {
    claim_level <- as.character(claim_level)
    if (identical(claim_level, "suggestive")) return("consistent_with")
    claim_level
}

# Attach the evidence a finding rests on, when the mapping is forced.
# See `.dev/spec-ai-evidence-authorship.md`.
#
# The problem this solves: `sclet_ai_derive_finding_claim_level()` is correct
# and the ceiling it consults is computed from evidence -- but with no
# `evidence_refs` on the finding there is nothing to reason about, so every
# finding lands at `hypothesis` through the "no refs" path. Measured on pbmc3k
# (`sclet-benchmark-ai/results/derived_claim_levels.jsonl`, 2026-10-13): the two
# `with_evidence` conditions had `derivation_available = TRUE`, the object
# carried `ev:findallmarkers` and the ceiling was `allowed`, and still 0 of 23
# and 0 of 13 findings moved -- because none cited anything.
#
# Why this is narrow on purpose: attaching evidence RAISES a claim, and an
# unsupported claim reported as clean is worse than the honest floor this
# replaces. The ceiling is the minimum over the cited nodes' own claim levels,
# so authorship cannot manufacture strength -- but it can still attach the
# WRONG node, and a finding about cell-cycle phase on an object whose only node
# is about markers would then look supported. So the rule attaches only when
# there is exactly one candidate and no choice to make, and otherwise leaves the
# finding alone. Under-attaching is the safe direction: a finding the engine
# cannot map keeps the floor, which is today's behaviour.
#
# Returns a list with the (possibly updated) `finding`, a logical `attached`,
# and the `rule` name that decided it. When `attached` is TRUE the finding's
# `evidence_refs` holds the sole node's id; the flag itself is not written onto
# the finding here -- `sclet_ai_recorded_findings()` and the audit entry are the
# places that record `evidence_refs_attached`.
sclet_ai_attach_finding_evidence <- function(object, finding) {
    if (!is.list(finding)) {
        return(list(finding = finding, attached = FALSE, rule = "none"))
    }
    # Two requested levels must not gain a reference, and the check has to come
    # before anything else:
    #
    #  - `causal` is banned outright. Attaching evidence to it would supply the
    #    very reference `validate_sclet_ai_result()` demands for a causal claim
    #    ("causal claim requires evidence_refs"), turning a result the engine
    #    rejects into one it accepts -- the ban defeated by helpfulness.
    #  - `estimated` is a provenance statement (a value was modelled, not
    #    measured), orthogonal to the evidence ladder. The derivation returns
    #    that verdict before consulting evidence at all, so a reference here
    #    would be attached to a finding whose level ignores it.
    requested <- as.character(finding$claim_level %||% "hypothesis")[[1L]]
    if (requested %in% c("causal", "estimated")) {
        return(list(finding = finding, attached = FALSE, rule = "requested_level_excluded"))
    }
    # The model's own citation is authoritative for what it cited. Never
    # overwrite it, even when a single node makes the engine's guess easy.
    existing <- as.character(
        if ("evidence_refs" %in% names(finding)) finding[["evidence_refs"]] else character()
    )
    existing <- existing[nzchar(existing)]
    if (length(existing)) {
        return(list(finding = finding, attached = FALSE, rule = "model_cited"))
    }
    if (!inherits(object, "SingleCellExperiment")) {
        return(list(finding = finding, attached = FALSE, rule = "none"))
    }
    nodes <- tryCatch(sclet_ai_evidence_get_all(object), error = function(e) list())
    if (!length(nodes)) {
        return(list(finding = finding, attached = FALSE, rule = "none"))
    }
    completed <- vapply(nodes, function(node) {
        identical(as.character(node$status %||% "completed"), "completed")
    }, logical(1L))
    nodes <- nodes[completed]
    # A user decision is not support for a finding. The ceiling already
    # separates these; excluding them here keeps the two rules aligned.
    is_decision <- vapply(nodes, function(node) {
        sclet_ai_evidence_kind_is_user_decision(node$kind %||% "")
    }, logical(1L))
    nodes <- nodes[!is_decision]
    if (length(nodes) != 1L) {
        # Zero candidates: nothing to attach. Two or more: attaching would be a
        # guess, and the honest answer is the floor.
        return(list(
            finding = finding,
            attached = FALSE,
            rule = if (!length(nodes)) "no_candidates" else "ambiguous_multiple_candidates"
        ))
    }
    finding$evidence_refs <- names(nodes)[[1L]]
    list(finding = finding, attached = TRUE, rule = "sole_evidence_node")
}

# Derive the claim level a finding should be recorded at, given the recorded
# evidence on the object. See `.dev/spec-ai-derived-claim-levels.md` for the
# derivation rule. The returned list has:
#   claim_level             - the level to record (or to refuse on)
#   requested_claim_level   - what the model asked for
#   derived_from_ceiling    - TRUE iff the ceiling decided the derivation
#   ceiling_status          - status of the underlying ceiling call (or
#                             "banned"/"provenance"/"no_refs" for the
#                             special-cased paths)
#   ceiling_reasons         - reasons from the ceiling function, so a level
#                             is never quoted without its basis
#   refused                 - TRUE iff this finding is a hard problem
#                             (currently only `causal`); the audit should
#                             never silently downgrade such a finding
sclet_ai_derive_finding_claim_level <- function(object, finding) {
    if (!is.list(finding)) {
        return(list(
            claim_level = "hypothesis",
            requested_claim_level = "hypothesis",
            derived_from_ceiling = FALSE,
            ceiling_status = "no_finding",
            ceiling_reasons = character(),
            refused = FALSE
        ))
    }
    requested <- as.character(finding$claim_level %||% "hypothesis")[[1L]]

    # 1. provenance: `estimated` is a statement about how a value was
    # produced, not a strength statement. It is orthogonal to the ceiling
    # and must survive the audit untouched.
    if (identical(requested, "estimated")) {
        return(list(
            claim_level = "estimated",
            requested_claim_level = "estimated",
            derived_from_ceiling = FALSE,
            ceiling_status = "provenance",
            ceiling_reasons = character(),
            refused = FALSE
        ))
    }

    # 2. causal: a hard ban. The audit must still report this as a problem
    # (`refused = TRUE`); the field-level derivation is `hypothesis` (the only
    # level a `causal` finding could ever land at) but the request is kept
    # visible via `requested_claim_level = "causal"`, and the audit entry
    # carries the causal problem. Per the spec, this is never a silent
    # downgrade: `RecordAIResult()` refuses on a `causal` finding because the
    # audit status is `overclaimed`.
    if (identical(requested, "causal")) {
        return(list(
            claim_level = "hypothesis",
            requested_claim_level = "causal",
            derived_from_ceiling = FALSE,
            ceiling_status = "banned",
            ceiling_reasons = "causal is not a claim level the AI may assert",
            refused = TRUE
        ))
    }

    # `[[` rather than `$`: a recorded finding also carries
    # `evidence_refs_attached` (logical), and `$` partial-matches, so a finding
    # with no refs would hand the ladder the flag `FALSE` as if it were an
    # evidence id -- which the ceiling then rejects as unknown, silently
    # pinning the finding to the floor for the wrong reason.
    refs <- if (is.list(finding) && "evidence_refs" %in% names(finding)) {
        finding[["evidence_refs"]]
    } else {
        NULL
    }
    refs <- as.character(refs %||% character())
    refs <- refs[nzchar(refs)]

    # 3a. no refs: the audit's previously-suggested floor was `observed`,
    # but the spec requires `hypothesis` as the correct floor for a finding
    # that asserts nothing without evidence. `hypothesis` is the only level
    # that explicitly self-limits and asserts nothing.
    if (!length(refs)) {
        return(list(
            claim_level = "hypothesis",
            requested_claim_level = requested,
            derived_from_ceiling = FALSE,
            ceiling_status = "no_refs",
            ceiling_reasons = "no evidence_refs supplied",
            refused = FALSE
        ))
    }

    # 3b. compute the ceiling against the cited evidence. An error (e.g., a
    # stale or unknown ref) is reported as `unsupported` and lands the
    # finding at `hypothesis` rather than crashing the audit.
    # `hypothesis` (and any other level outside the evidence ladder) is not a
    # valid `proposed_claim_level` for the ceiling, which validates against
    # `sclet_ai_claim_levels`. Passing it anyway made the ceiling raise
    # "proposed_claim_level must be one of: ...", which the tryCatch below
    # swallowed as `unsupported` and landed the finding back at `hypothesis`.
    # Since the model emits `hypothesis` for every finding, that silent
    # fallback made the whole derivation a no-op on real data. A requested
    # level outside the ladder asserts nothing (rank 0), so passing NULL lets
    # the ceiling simply report what the evidence can carry.
    proposed <- sclet_ai_finding_allowed_level(requested)
    if (is.na(sclet_ai_claim_rank(proposed))) {
        proposed <- NULL
    }
    ceiling <- tryCatch(
        sclet_ai_claim_ceiling(
            object,
            evidence_ids = refs,
            proposed_claim_level = proposed
        ),
        error = function(e) list(
            status = "unsupported",
            ceiling_claim_level = NULL,
            n_independent_groups = 0L,
            reasons = conditionMessage(e)
        )
    )

    if (identical(ceiling$status, "unsupported")) {
        return(list(
            claim_level = "hypothesis",
            requested_claim_level = requested,
            derived_from_ceiling = FALSE,
            ceiling_status = "unsupported",
            ceiling_reasons = as.character(ceiling$reasons %||% "no independent support"),
            refused = FALSE
        ))
    }

    # 3c. The derived level is a function of the EVIDENCE, not a copy of the
    # model's word. The ceiling is the strongest level the cited evidence can
    # carry, and it is computed as the minimum over the supporting nodes' own
    # claim levels -- so taking it can never assert more than the evidence
    # already declares for itself.
    #
    # It is deliberately NOT min(requested, ceiling): capping at the request
    # would let a model that always answers `hypothesis` pin every finding to
    # the floor, which is exactly the state this change exists to fix (114 of
    # 114 measured findings were `hypothesis`). A request BELOW the ceiling is
    # raised to the ceiling; a request ABOVE it is lowered to the ceiling. The
    # request is always preserved in `requested_claim_level`, so the movement
    # stays visible either way.
    ceiling_rank <- sclet_ai_claim_rank(ceiling$ceiling_claim_level)
    if (is.na(ceiling_rank)) ceiling_rank <- 0L
    derived <- if (ceiling_rank < 1L) "hypothesis" else sclet_ai_claim_from_rank(ceiling_rank)
    list(
        claim_level = derived,
        requested_claim_level = requested,
        derived_from_ceiling = TRUE,
        ceiling_status = ceiling$status,
        ceiling_reasons = as.character(ceiling$reasons %||% character()),
        refused = FALSE
    )
}

#' Audit the findings of an AI result against the recorded evidence
#'
#' `validate_sclet_ai_result()` can only check a result's shape, because it has
#' no object to look at. This function performs the object-aware half of the
#' contract: for every finding, it derives the claim level the finding should
#' be recorded at from the recorded evidence, and reports the per-finding
#' verdict.
#'
#' Three rules are enforced. First, `"causal"` is never an acceptable AI
#' claim level: the advanced-analysis contract reserves causal interpretation
#' for a separately designed and separately evidenced causal study, and no
#' amount of observational evidence lifts an AI finding into it. A `causal`
#' finding is reported as `overclaimed` and is never silently downgraded.
#' Second, a finding whose requested level exceeds the evidence ceiling is
#' reported as `downgraded` with the corrected level attached, so the
#' downgrade is visible rather than silent. Third, a finding that cites no
#' evidence (or unsupported evidence) is also reported as `downgraded`, with
#' the derived level set to `hypothesis` -- the only level that explicitly
#' self-limits and asserts nothing.
#'
#' The function only audits. It never rewrites the caller's `result`, never
#' downgrades a `causal` finding into silence, and never records anything.
#' `RecordAIResult()` is the function that applies the derived level to the
#' recorded finding.
#'
#' @param object A `SingleCellExperiment` object.
#' @param result An `sclet_ai_result` or a list with a `findings` element.
#' @return A list with an overall `status` of `"ok"`, `"downgraded"` or
#'   `"overclaimed"`, the per-finding `findings` audit, the `n_problems` and
#'   `n_downgraded` counts, and a `problems` character vector. Each
#'   `findings[[i]]` carries `claim_level` (derived), `requested_claim_level`
#'   (what the model asked), `derived_from_ceiling` (logical), `ceiling_reasons`
#'   (the basis), `allowed_claim_level` (kept for backward compatibility;
#'   equal to the derived `claim_level`), and a per-entry `status`.
#' @noRd
sclet_ai_audit_result_claims <- function(object, result) {
    findings <- if (is.list(result)) result$findings else NULL
    if (!is.list(findings) || !length(findings)) {
        return(list(
            status = "ok", findings = list(),
            n_problems = 0L, n_downgraded = 0L, problems = character()
        ))
    }
    problems <- character()
    audit <- list()
    n_downgraded <- 0L
    for (i in seq_along(findings)) {
        finding <- findings[[i]]
        requested <- if (is.list(finding)) as.character(finding$claim_level %||% "hypothesis")[[1L]] else "hypothesis"
        # Authorship runs first: a finding that gains a reference gains its
        # derived level in the same pass, because the derivation reads refs.
        attachment <- sclet_ai_attach_finding_evidence(object, finding)
        finding <- attachment$finding
        derivation <- sclet_ai_derive_finding_claim_level(object, finding)
        derived <- derivation$claim_level

        # Per-entry verdict:
        #  - `causal` -> overclaimed (hard ban)
        #  - derivation differs from request -> downgraded (correction)
        #  - otherwise -> ok
        entry_status <- "ok"
        if (isTRUE(derivation$refused)) {
            entry_status <- "overclaimed"
            problems <- c(problems, paste0(
                "findings[[", i, "]]: 'causal' is not a claim level the AI may assert; ",
                "causal interpretation requires a separately designed and evidenced causal study. ",
                "Use 'hypothesis' at most."
            ))
        } else if (!identical(requested, derived)) {
            entry_status <- "downgraded"
            n_downgraded <- n_downgraded + 1L
        }

        entry <- list(
            index = i,
            claim_level = derived,
            requested_claim_level = requested,
            derived_from_ceiling = isTRUE(derivation$derived_from_ceiling),
            allowed_claim_level = derived,
            ceiling_status = derivation$ceiling_status,
            ceiling_reasons = derivation$ceiling_reasons,
            # Authorship, made visible rather than laundered into the verdict.
            # A reader must be able to tell a reference the model cited from one
            # the engine attached on its behalf, because only the former is
            # evidence the model actually reasoned about.
            evidence_refs_attached = isTRUE(attachment$attached),
            attachment_rule = as.character(attachment$rule %||% "none"),
            status = entry_status,
            problems = character()
        )
        audit[[length(audit) + 1L]] <- entry
    }
    overall_status <- if (length(problems)) "overclaimed"
        else if (n_downgraded > 0L) "downgraded"
        else "ok"
    list(
        status = overall_status,
        findings = audit,
        n_problems = length(problems),
        n_downgraded = n_downgraded,
        problems = problems
    )
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