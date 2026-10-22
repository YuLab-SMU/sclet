## Shared AI-layer fixtures.
##
## Loaded by testthat for every test file (and for a single test_file() run),
## so fixture definitions live here exactly once instead of being duplicated
## or rediscovered by scanning another test file's source.

## The per-cell-vector guard: TRUE when any atomic leaf of `x` holds one entry
## per cell, which is the shape a raw per-cell vector would take if it leaked
## into a bounded evidence node. Recursion covers unnamed list elements too, so
## a structure that buries a vector one level deeper is still caught.
sclet_ai_test_has_per_cell <- function(x, n_cells) {
    if (is.list(x)) {
        any(vapply(x, sclet_ai_test_has_per_cell, logical(1L), n_cells = n_cells))
    } else {
        is.atomic(x) && length(x) >= n_cells
    }
}

## ---------------------------------------------------------------------------
## Cell-cell interaction fixtures
## ---------------------------------------------------------------------------

sclet_ai_test_cci_genes <- function() {
    c("CXCL12", "ACKR1", "HLA-A", "CD8A", "APP", "CD74", "ACTB", "GAPDH")
}

sclet_ai_test_cci_object <- function(n_cells = 40L) {
    genes <- sclet_ai_test_cci_genes()
    counts <- matrix(rpois(length(genes) * n_cells, lambda = 8L),
        nrow = length(genes), ncol = n_cells,
        dimnames = list(genes, paste0("cell_", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    ## Labels that are NOT already anonymization-shaped, so a leak of a raw
    ## label and a legitimate code are distinguishable in the assertions below.
    sce$cluster <- factor(rep(c("beta cells", "alpha cells"), length.out = n_cells),
        levels = c("beta cells", "alpha cells"))
    ActiveIdent(sce) <- "cluster"
    sce
}

## The standardized interaction table the CCI runners persist. Columns and
## their order mirror R/cci.R (source_group, target_group, ligand, receptor,
## pathway, probability, p_value, method).
sclet_ai_test_cci_table <- function(method = "CellChat") {
    data.frame(
        source_group = c("beta cells", "beta cells", "alpha cells", "alpha cells", "beta cells"),
        target_group = c("alpha cells", "alpha cells", "beta cells", "beta cells", "beta cells"),
        ligand = c("CXCL12", "APP", "HLA-A", "CXCL12", "ACKR1"),
        receptor = c("ACKR1", "CD74", "CD8A", "ACKR1", "CXCL12"),
        pathway = c("CXCL12_ACKR1", "APP_CD74", "HLA-A_CD8A", "CXCL12_ACKR1", "ACKR1_CXCL12"),
        probability = c(0.91, 0.62, 0.44, 0.31, 0.08),
        p_value = c(0.001, 0.02, 0.34, 0.48, 0.99),
        method = method,
        stringsAsFactors = FALSE
    )
}

## Persist a communication result the way the real runners do: the interaction
## table into the analysis layer (where get_cci() reads it from) and a
## type="communication" state record that carries only counts.
sclet_ai_test_cci_store <- function(object, id = "cci1", method = "CellChat",
                                   analysis_key = "cellchat", table = sclet_ai_test_cci_table(method),
                                   group = "cluster") {
    object <- sclet_set_analysis(object, analysis_key, list(
        id = id,
        method = method,
        object = NULL,
        interactions = table,
        group = group
    ))
    sclet_set_analysis_state(
        object = object,
        type = "communication",
        id = id,
        method = method,
        inputs = list(assay = "counts", layer = NULL, group = group),
        artifacts = list(analysis_key = analysis_key),
        params = list(species = "human"),
        summary = list(n_groups = length(unique(table$source_group))),
        active = TRUE
    )
}

## Stub for RunCCI(). It writes through the real storage path so everything
## downstream is exercised; `writes = FALSE` reproduces the silent no-record
## case the action has to detect.
sclet_ai_test_cci_engine_stub <- function(id = "cci1", method = "CellChat",
                                          analysis_key = "cellchat", writes = TRUE,
                                          table = sclet_ai_test_cci_table(method)) {
    function(object, method = NULL, ...) {
        if (!writes) return(object)
        sclet_ai_test_cci_store(object, id = id, method = method,
            analysis_key = analysis_key, table = table)
    }
}

## Pin the CCI engine gate for the duration of the calling test. The gate asks
## the real library, so an unpinned test asserts whatever the machine happens to
## have installed: CI installs nichenetr through its own setup step, while many
## developer machines do not have it, and CellChat is absent from some. The gate
## is what is under test here, not this environment, so both branches are forced
## deterministically: CellChat and basilisk present, nichenetr absent. Every
## other package keeps its real answer, and the pin is undone when the test ends.
sclet_ai_test_cci_pin_engines <- function(env = parent.frame()) {
    real_req_ns <- base::requireNamespace
    testthat::local_mocked_bindings(
        requireNamespace = function(package, ...) {
            if (identical(package, "CellChat")) return(TRUE)
            if (identical(package, "basilisk")) return(TRUE)
            if (identical(package, "nichenetr")) return(FALSE)
            real_req_ns(package, ...)
        },
        .package = "base",
        .env = env
    )
}

sclet_ai_test_cci_registry <- function(object) {
    AIDefaultExecutionRegistry(object, include = c("read", "communication"))
}

sclet_ai_test_cci_plan <- function(object, registry, params) {
    plan <- sclet:::new_sclet_ai_plan(
        task = "communication",
        context_fingerprint = GetAnalysisLedger(object)$fingerprint,
        actions = list(list(id = "cc1", action = "run_cci", params = as.list(params)))
    )
    ValidateAIPlan(plan, object = object, registry = registry, strict = FALSE)
}

## Every character value anywhere in a nested object.
sclet_ai_test_cci_characters <- function(x) {
    if (is.character(x)) return(x)
    if (is.list(x)) return(unlist(lapply(x, sclet_ai_test_cci_characters), use.names = FALSE))
    character()
}

## Every piece of text an object can carry, as ONE string (a character vector of
## length 1). capture.output(str()) is NOT used for the leak scan: str() elides
## long character vectors and truncates wide data frames, so a value that really
## leaked could sit in an omitted section and the scan would pass vacuously.
## Names are included because a field called `ligand` is itself a leak signal.
##
## The return contract is length 1 on purpose: every caller feeds the result to
## grepl(..., fixed = TRUE), and a per-leaf character vector would make that
## grepl return one logical per leaf instead of one answer for the whole scan.
##
## Pure recursion, no closure or `<<-`: this is a shared helper now, and an
## accumulator held in an enclosing environment is the kind of state a later
## caller re-uses by mistake.
sclet_ai_test_cci_all_text <- function(x) {
    pieces <- character()
    if (is.data.frame(x)) {
        for (column in colnames(x)) {
            pieces <- c(pieces, column)
            pieces <- c(pieces, sclet_ai_test_cci_all_text(x[[column]]))
        }
    } else if (is.list(x)) {
        if (!is.null(names(x))) pieces <- c(pieces, names(x))
        for (element in x) pieces <- c(pieces, sclet_ai_test_cci_all_text(element))
    } else if (is.character(x)) {
        pieces <- c(pieces, paste(x, collapse = "\r"))
    }
    paste(pieces, collapse = "\n")
}

## The identity strings this layer must never let escape.
sclet_ai_test_cci_identity_strings <- function(table) {
    unique(c(
        "CXCL12", "ACKR1", "HLA-A", "CD8A", "APP", "CD74",
        table$ligand, table$receptor, table$pathway,
        "beta cells", "alpha cells"
    ))
}
