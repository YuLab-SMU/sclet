## Evidence-node census driver (offline, deterministic pipeline only).
## Usage: Rscript .dev/census-run.R <stage>
##   stage: core|hvvg|scale|pca|graph|clust|umap|de|ann|dbl|rare|traj|vel|final
## Each stage runs EXACTLY ONE action through ExecuteAIPlan() and dumps the node
## set before and after. State: /tmp/sclet_census/sce.rds

suppressWarnings(suppressMessages({
    library(pkgload)
    load_all("/home/wang/data/source/omics/sclet", export_all = TRUE, quiet = TRUE)
}))

STATEDIR <- "/tmp/sclet_census"
dir.create(STATEDIR, showWarnings = FALSE, recursive = TRUE)
RDS <- file.path(STATEDIR, "sce.rds")

stage <- commandArgs(trailingOnly = TRUE)[[1]]

cat_value <- function(x) {
    if (is.null(x)) return("NULL")
    paste(format(x), collapse = ", ")
}

node_ids <- function(object) {
    tryCatch(names(sclet:::sclet_ai_evidence_get_all(object)), error = function(e) character())
}

dump_nodes <- function(object, label) {
    nodes <- sclet:::sclet_ai_evidence_get_all(object)
    ids <- names(nodes)
    lines <- paste0("=== DUMP ", label, " | n_nodes=", length(ids))
    if (!length(ids)) {
        lines <- c(lines, "(no evidence nodes)")
    }
    for (i in seq_along(ids)) {
        id <- ids[[i]]; n <- nodes[[i]]
        lines <- c(lines, paste0("--- node ", i, " id=", id))
        for (f in c("kind","claim_level","source","status","dependency_group",
                    "object_fingerprint","created_at")) {
            lines <- c(lines, paste0("    ", f, " = ", cat_value(n[[f]])))
        }
        lines <- c(lines, paste0("    parents = ", cat_value(n$parents)))
        sc <- n$scope
        if (is.list(sc) && length(sc)) {
            lines <- c(lines, paste0("    scope = list(",
                paste(sprintf("%s=%s", names(sc), vapply(sc, cat_value, character(1))), collapse = ", "), ")"))
        } else {
            lines <- c(lines, paste0("    scope = ", cat_value(sc)))
        }
        v <- n$values
        if (is.list(v) && length(v)) {
            summ <- vapply(names(v), function(nm) {
                el <- v[[nm]]
                if (is.list(el)) sprintf("%s=list(len=%d,names=%s)", nm, length(el), paste(names(el) %||% rep("", length(el)), collapse = ","))
                else sprintf("%s=%s(len=%d,first=%s)", nm, class(el)[[1]], length(el), cat_value(utils::head(el, 2)))
            }, character(1))
            lines <- c(lines, paste0("    values names(", length(v), "): ", paste(names(v), collapse = ", ")))
            lines <- c(lines, paste0("    values shape: ", paste(summ, collapse = " | ")))
        } else {
            lines <- c(lines, paste0("    values = ", cat_value(v)))
        }
        lines <- c(lines, paste0("    all fields: ", paste(names(n), collapse = ", ")))
    }
    out <- paste(lines, collapse = "\n")
    writeLines(out, file.path(STATEDIR, paste0("dump-", label, ".txt")))
    cat(out, "\n", sep = "")
    invisible(ids)
}

ceiling_report <- function(object, label) {
    cl <- tryCatch(sclet:::sclet_ai_claim_ceiling(object), error = function(e) e)
    if (inherits(cl, "error")) {
        cat("=== CEILING ", label, " ERROR: ", conditionMessage(cl), "\n", sep = "")
        return(invisible(NULL))
    }
    txt <- c(paste0("=== CEILING ", label),
        paste0("    status = ", cat_value(cl$status)),
        paste0("    ceiling_claim_level = ", cat_value(cl$ceiling_claim_level)),
        paste0("    n_nodes = ", cat_value(cl$n_nodes)),
        paste0("    n_independent_groups = ", cat_value(cl$n_independent_groups)),
        paste0("    reasons = ", paste(unlist(cl$reasons), collapse = " || ")),
        paste0("    n_dependent_pairs = ", length(cl$dependent_pairs)),
        paste0("    groups: ", paste(vapply(cl$groups %||% list(), function(g)
            sprintf("#%d n=%d nodes=[%s] sources=[%s]", g$group, g$n_nodes,
                    paste(g$nodes, collapse = ","), paste(g$sources, collapse = ",")),
            character(1)), collapse = " ;; ")),
        paste0("    user_decision_nodes = ", cat_value(cl$user_decision_nodes)))
    out <- paste(txt, collapse = "\n")
    writeLines(out, file.path(STATEDIR, paste0("ceiling-", label, ".txt")))
    cat(out, "\n", sep = "")
    invisible(NULL)
}

## One action, one plan, one ExecuteAIPlan call.
run_one <- function(sce, action, params = list(), label, retries = 0L) {
    before <- node_ids(sce)
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    if (!action %in% names(reg)) {
        cat("### ", label, ": action '", action, "' NOT IN REGISTRY\n", sep = "")
        return(sce)
    }
    plan <- new_sclet_ai_plan(task = paste0("census_", label),
        actions = list(list(id = paste0("s_", label), action = action,
                            params = params, requires_confirmation = TRUE,
                            max_retries = retries)),
        requires_confirmation = TRUE)
    val <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    cat("--- ", label, ": validation valid=", isTRUE(val$valid), "\n", sep = "")
    if (length(val$errors))   cat("      errors:   ", paste(val$errors, collapse = " | "), "\n", sep = "")
    if (length(val$warnings)) cat("      warnings: ", paste(val$warnings, collapse = " | "), "\n", sep = "")
    if (!isTRUE(val$valid)) {
        cat("### ", label, ": NOT EXECUTED (plan invalid; errors quoted verbatim above)\n", sep = "")
        return(sce)
    }
    t0 <- Sys.time()
    ex <- tryCatch(ExecuteAIPlan(object = sce, plan = plan, registry = reg,
        validation = val, dry_run = FALSE,
        confirmation = val$confirmation_token, record = FALSE),
        error = function(e) e)
    if (inherits(ex, "error")) {
        cat("### ", label, ": ExecuteAIPlan raised: ", conditionMessage(ex), "\n", sep = "")
        return(sce)
    }
    for (r in ex$results) {
        cat(sprintf("      step %s: status=%s dur=%.1fs\n", r$id, r$status,
                    as.numeric(r$duration_sec %||% NA)))
        if (!is.null(r$error)) cat("      step error: ", r$error, "\n", sep = "")
    }
    cat("      overall status = ", ex$status, " (elapsed ",
        format(difftime(Sys.time(), t0, units = "secs")), ")\n", sep = "")
    out <- ex$object
    after <- node_ids(out)
    new <- setdiff(after, before)
    cat(sprintf("### NODE COUNT %s: before=%d after=%d new=%d [%s]\n",
                label, length(before), length(after), length(new),
                if (length(new)) paste(new, collapse = ", ") else "-"))
    if (length(new)) dump_nodes(out, paste0(label, "_new"))
    ceiling_report(out, label)
    out
}

if (!file.exists(RDS)) stop("no state at ", RDS)
sce <- readRDS(RDS)
cat("state loaded: dim", paste(dim(sce), collapse = " x "),
    "| assays", paste(assayNames(sce), collapse = ","),
    "| nodes", length(node_ids(sce)), "\n")

sce <- switch(stage,
    core  = run_one(sce, "normalize_data", list(), "core"),
    hvg   = run_one(sce, "find_variable_features", list(), "hvg"),
    scale = run_one(sce, "scale_data", list(), "scale"),
    pca   = run_one(sce, "run_pca", list(ncomponents = 50L), "pca"),
    graph = run_one(sce, "find_neighbors", list(dims = as.integer(1:20), k = 10L), "graph"),
    clust = run_one(sce, "find_clusters", list(resolution = 0.5), "clust"),
    umap  = run_one(sce, "run_umap", list(), "umap"),
    de    = run_one(sce, "run_de_test", list(), "de"),
    ann   = run_one(sce, "run_annotation",
                    list(ref = "HumanPrimaryCellAtlasData", labels = "label.cls",
                         name = "hpcapca"), "ann"),
    dbl   = run_one(sce, "run_doublet_detection", list(), "dbl"),
    rare  = run_one(sce, "run_rare_cell_detection",
                    list(reduction = "PCA", dims = as.list(1:20),
                         rare_threshold = 100L, name = "rareq"), "rare"),
    traj  = {
        labs <- sort(unique(as.character(Idents(sce))))
        cols <- colnames(colData(sce))
        gcol <- intersect(c("clusters", "cluster", "ident", "seurat_clusters", "labels"), cols)[[1]]
        if (is.na(gcol)) gcol <- "ident"
        cat("      traj params: group=", gcol, " start_cluster=", labs[[1]],
            " labels=", paste(labs, collapse = ","), "\n", sep = "")
        run_one(sce, "run_trajectory",
                list(group = gcol, start_cluster = labs[[1]], reduction = "UMAP",
                     name = "slingshot", seed = 2025L), "traj")
    },
    vel   = run_one(sce, "run_velocity", list(mode = "deterministic", name = "velocity"), "vel"),
    final = { dump_nodes(sce, "99_final"); ceiling_report(sce, "99_final")
              cat("colData:", paste(colnames(colData(sce)), collapse = ","), "\n")
              cat("reducedDims:", paste(reducedDimNames(sce), collapse = ","), "\n")
              cat("assays:", paste(assayNames(sce), collapse = ","), "\n")
              cat("graphs:", paste(names(sclet_get_state(sce)$graphs %||% list()), collapse = ","), "\n")
              sce },
    stop("unknown stage: ", stage))

saveRDS(sce, RDS)
cat("[state saved; nodes now = ", length(node_ids(sce)), "]\n", sep = "")
