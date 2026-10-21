`%||%` <- function(x, y) if (is.null(x)) y else x
# Build a per-function example plan for the Rd files that have no \examples.
#
# For every exported function it reports:
#   * the R source file that defines it
#   * the Suggests packages the function body may need (via requireNamespace)
#   * whether the body reads recorded state (i.e. an example needs a pipeline
#     or a pre-built object before the call is meaningful)
#
# Output: a TSV with one row per target Rd, consumed by the example-writing
# agents. This is a hint file, not an authority -- the agent still has to read
# the function.

suppressWarnings(suppressMessages(library(methods)))

wd <- getwd()
desc <- readLines(file.path(wd, "DESCRIPTION"), warn = FALSE)
field <- function(name) {
    i <- grep(paste0("^", name, ":"), desc)
    if (!length(i)) return(character())
    start <- i[[1]] + 1L
    if (start > length(desc)) return(character())
    tail <- desc[start:length(desc)]
    stop_at <- grep("^[A-Za-z][A-Za-z0-9_.]*:", tail)
    if (length(stop_at)) tail <- tail[seq_len(stop_at[[1]] - 1L)]
    # One package per line; a trailing comma belongs to the line, not the
    # next one. Split on newlines FIRST -- collapsing ",\n" into " " would
    # merge the whole block into a single string.
    parts <- trimws(tail)
    parts <- sub(",$", "", parts)
    parts <- parts[nzchar(parts)]
    parts
}
suggests <- field("Suggests")
imports <- field("Imports")
cat("Suggests:", length(suggests), " Imports:", length(imports), "\n")

ns <- readLines(file.path(wd, "NAMESPACE"), warn = FALSE)
exports <- gsub("^export\\(|\\)$", "", grep("^export\\(", ns, value = TRUE))
exports <- gsub('"', "", exports)
exports <- sort(unique(exports))

# Parse every R file once and keep the source of each top-level assignment.
src_files <- list.files(file.path(wd, "R"), pattern = "[.]R$", full.names = TRUE)
bodies <- list()
for (f in src_files) {
    exprs <- tryCatch(parse(f, keep.source = FALSE), error = function(e) NULL)
    if (is.null(exprs)) next
    for (e in exprs) {
        if (is.call(e) && length(e) >= 3L &&
            identical(as.character(e[[1]]), "<-") && is.name(e[[2]])) {
            nm <- as.character(e[[2]])
            txt <- paste(deparse(e), collapse = "\n")
            if (!is.null(bodies[[nm]])) next
            bodies[[nm]] <- list(file = basename(f), text = txt)
        }
    }
}
cat("parsed assignments:", length(bodies), "\n")

# Does this body mention a Suggests package (by requireNamespace or pkg::)?
suggests_needed <- function(txt) {
    hits <- character()
    for (p in suggests) {
        if (grepl(paste0("requireNamespace\\(\\s*[\"']", p, "[\"']"), txt)) hits <- c(hits, p)
    }
    hits <- unique(hits)
    if (!length(hits)) {
        for (p in suggests) {
            if (grepl(paste0("(?<![A-Za-z0-9._])", p, "::"), txt, perl = TRUE)) hits <- c(hits, p)
        }
        hits <- unique(hits)
    }
    hits
}

# Does the body read recorded state / a previous analysis result?
state_readers <- c(
    "get_detest", "get_enrichment", "get_milo", "get_cci", "get_cellchat",
    "get_cellphonedb", "get_cellrank", "get_annotation", "get_mapping",
    "get_hvg", "get_graph", "get_batch", "get_geneset_scoring",
    "get_colocalization", "sclet_get_state_record", "sclet_get_analysis",
    "sclet_get_active_state", "sclet_get_active_ident",
    "ActiveIdent", "DefaultReduction", "DefaultAssay"
)
reads_state <- function(txt) {
    pats <- vapply(state_readers, function(f) {
        grepl(paste0("(?<![A-Za-z0-9._])", f, "\\s*\\("), txt, perl = TRUE)
    }, logical(1))
    any(pats)
}
# Does the body need a python / basilisk environment?
needs_python <- function(txt) {
    grepl("sclet_basilisk_run|sclet_env_run|sclet_env_path", txt)
}

target_rd <- commandArgs(trailingOnly = TRUE)
if (!length(target_rd)) {
    target_rd <- list.files(file.path(wd, "man"), pattern = "[.]Rd$", full.names = TRUE)
}
target_rd <- target_rd[!grepl("sclet_ai_", target_rd)]

# An Rd file name is NOT the object name. `subset-cell.Rd` documents
# `subset_cell`; `reexports.Rd` documents a dozen re-exported names; a
# `@rdname` topic page documents several functions at once. So resolve the
# real object names from the Rd's own \name / \alias tags, and only then look
# up which R/ file defines them.
rd_tags <- function(rd, tag) {
    rdo <- tryCatch(tools::parse_Rd(rd), error = function(e) NULL)
    if (is.null(rdo)) return(character())
    out <- character()
    for (part in rdo) {
        if (!identical(attr(part, "Rd_tag"), paste0("\\", tag))) next
        out <- c(out, paste(as.character(part), collapse = ""))
    }
    out <- trimws(out)
    out[nzchar(out)]
}
# Topic pages that document several objects at once, or objects that are not
# plain assignments (the re-exports are `pkg::name` expressions).
topic_owner <- c(
    "reexports.Rd" = "reexport.R",
    "sclet-package.Rd" = "sclet-package.R",
    "advanced-ai-diagnostics.Rd" = "ai-diagnostics.R"
)

out <- file.path(wd, ".dev", "tools", "example-plan.tsv")
con <- file(out, "w")
writeLines("rd\tname\tsrc_file\tsuggests_needed\treads_state\tneeds_python", con)
n <- 0L
for (rd in target_rd) {
    txt <- paste(readLines(rd, warn = FALSE), collapse = "\n")
    if (grepl("\\\\examples", txt)) next
    base <- basename(rd)
    nm <- sub("[.]Rd$", "", base)
    # candidate object names for this page
    cands <- unique(c(rd_tags(rd, "name"), rd_tags(rd, "alias")))
    cands <- cands[nzchar(cands)]
    b <- NULL
    for (cn in cands) {
        if (!is.null(bodies[[cn]])) {
            b <- bodies[[cn]]
            nm <- cn
            break
        }
    }
    if (is.null(b) && base %in% names(topic_owner)) {
        b <- list(file = topic_owner[[base]], text = "")
        nm <- cands[[1]] %||% nm
    }
    if (is.null(b)) {
        writeLines(paste(base, nm, "UNKNOWN", "", "", "", sep = "\t"), con)
        n <- n + 1L
        next
    }
    sg <- suggests_needed(b$text)
    writeLines(
        paste(base, nm, b$file,
              if (length(sg)) paste(sg, collapse = ",") else "-",
              if (reads_state(b$text)) "yes" else "no",
              if (needs_python(b$text)) "yes" else "no",
              sep = "\t"),
        con
    )
    n <- n + 1L
}
close(con)
cat("wrote", out, "with", n, "rows\n")

