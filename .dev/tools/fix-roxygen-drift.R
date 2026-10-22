# Repair the roxygen/ Rd drift: nine functions whose roxygen block in R/ is
# older than the man/ page that was committed.
#
# The committed Rd has richer @param text than the roxygen source, so running
# roxygenise() as-is DELETES that documentation -- and then R CMD check's
# codoc starts reporting undocumented arguments. This script copies the
# missing parameter descriptions out of the committed Rd and into the roxygen
# block, so the next roxygenise() reproduces the committed text.
#
# The committed Rd is read from git HEAD, not from the working tree, so the
# script is idempotent and does not depend on a roxygenise having run.
#
#   Rscript .dev/tools/fix-roxygen-drift.R            # dry run
#   Rscript .dev/tools/fix-roxygen-drift.R --apply
#
# Two things this script got wrong the first time, both recorded here because
# the failure mode is silent:
#
#   * The Rd's text was read with `as.character()` on the parsed fragment. On
#     an Rd fragment that DEPARSES: `\code{colnames(expr)}` comes back as the
#     string `list("colnames(expr)")`. The repaired blocks were filled with
#     that junk. The text is now extracted structurally, tag by tag.
#   * Parameters were matched by their first name. Roxygen allows a combined
#     `@param min_expr,min_frac` tag, and the Rd it produces has ONE
#     \item{min_expr, min_frac}. Matching on the first name alone made an
#     already-documented parameter look missing, and the script ADDED a
#     second one -- which R CMD check then reported as a duplicated
#     \argument. Parameters are now matched on the SET of names a tag covers.

wd <- getwd()
apply <- "--apply" %in% commandArgs(trailingOnly = TRUE)

# ---------------------------------------------------------------- Rd -> text

# An Rd fragment to plain text, tag by tag. `as.character()` on a fragment
# deparses it, which is not what a human wrote.
rd_fragment_text <- function(x) {
    if (is.character(x)) return(paste(x, collapse = ""))
    if (!is.list(x)) return("")
    tag <- attr(x, "Rd_tag")
    inner <- paste(vapply(x, rd_fragment_text, character(1)), collapse = "")
    if (is.null(tag) || length(tag) != 1L) return(inner)
    switch(tag,
        "\\code"   = paste0("`", inner, "`"),
        "\\link"   = paste0("[", inner, "]"),
        "\\emph"   = paste0("*", inner, "*"),
        "\\strong" = paste0("**", inner, "**"),
        "\\pkg"    = paste0("`", inner, "`"),
        "\\dots"   = "...",
        inner)
}

# The \arguments{\item{names}{text} ...} entries of an Rd file, keyed by the
# SET of names each entry covers.
rd_args <- function(path) {
    rdo <- tryCatch(tools::parse_Rd(path), error = function(e) NULL)
    if (is.null(rdo)) return(list())
    out <- list()
    for (part in rdo) {
        if (!identical(attr(part, "Rd_tag"), "\\arguments")) next
        for (item in as.list(part)) {
            if (!identical(attr(item, "Rd_tag"), "\\item")) next
            kids <- as.list(item)
            nm <- paste(as.character(kids[[1]]), collapse = "")
            txt <- rd_fragment_text(kids[[2]])
            out[[param_key(nm)]] <- list(names = param_names(nm), text = trimws(txt))
        }
    }
    out
}

# Normalise a parameter-name list to a comparable key: names are split on
# commas, trimmed and sorted, so `@param min_expr,min_frac` and
# `\item{min_expr, min_frac}` land on the same key.
param_names <- function(nm) {
    nms <- trimws(strsplit(as.character(nm), ",")[[1]])
    nms <- nms[nzchar(nms)]
    sort(unique(nms))
}
param_key <- function(nm) paste(param_names(nm), collapse = ",")

# ---------------------------------------------------------- roxygen blocks

# The @param entries of a roxygen block, keyed the same way.
block_params <- function(block) {
    out <- list()
    i <- 1
    while (i <= length(block)) {
        m <- regmatches(block[i], regexec("^\\s*#'\\s*@param\\s+([^[:space:]]+)(.*)$", block[i]))[[1]]
        if (!length(m)) {
            i <- i + 1
            next
        }
        names_raw <- m[[2]]
        rest <- trimws(m[[3]])
        text <- rest
        j <- i + 1
        while (j <= length(block) && grepl("^\\s*#'", block[[j]]) &&
               !grepl("^\\s*#'\\s*@", block[[j]])) {
            text <- paste(text, trimws(sub("^\\s*#'\\s?", "", block[[j]])))
            j <- j + 1L
        }
        key <- param_key(names_raw)
        out[[key]] <- list(names = param_names(names_raw), text = trimws(text),
                           first_line = i, last_line = j - 1L, raw_names = names_raw)
        i <- j
    }
    out
}

# Write one @param entry, wrapped the way the package's blocks do.
wrap_param <- function(nm, txt) {
    prefix <- paste0("#' @param ", nm, " ")
    width <- 78L
    words <- strsplit(paste(txt, collapse = " "), "[[:space:]]+")[[1]]
    if (!length(words)) return(paste0(prefix, txt))
    lines <- character()
    cur <- prefix
    for (w in words) {
        if (nchar(cur) + nchar(w) + 1L > width && nchar(cur) > nchar(prefix)) {
            lines <- c(lines, cur)
            cur <- paste0("#'    ", w)
        } else {
            cur <- paste0(cur, if (identical(cur, prefix)) "" else " ", w)
        }
    }
    c(lines, cur)
}

# ------------------------------------------------------------------- targets

targets <- list(
    list(rd = "man/CallMalignantCells.Rd",     src = "R/malignancy-combine.R",  fn = "CallMalignantCells"),
    list(rd = "man/CallPhenotypeCells.Rd",      src = "R/phenotype-evidence.R",  fn = "CallPhenotypeCells"),
    list(rd = "man/FindClusters.Rd",            src = "R/clustering.R",          fn = "FindClusters"),
    list(rd = "man/FindMalignantSeedGenes.Rd",  src = "R/malignancy.R",          fn = "FindMalignantSeedGenes"),
    list(rd = "man/GetAnalysisLedger.Rd",       src = "R/ai-context.R",          fn = "GetAnalysisLedger"),
    list(rd = "man/PhenotypeEvidenceProviders.Rd", src = "R/phenotype-evidence.R", fn = "PhenotypeEvidenceProviders"),
    list(rd = "man/RunCNVInference.Rd",         src = "R/malignancy-cnv.R",      fn = "RunCNVInference"),
    list(rd = "man/RunDEtest.Rd",               src = "R/markers.R",             fn = "RunDEtest"),
    list(rd = "man/sclet-package.Rd",           src = "R/sclet-package.R",       fn = NULL)
)

# A description is "materially thinner" when the committed Rd says at least
# this many more characters. Measured on CallMalignantCells: a mere markup
# difference (a bare `SingleCellExperiment` vs the same word in \code{})
# differs by 6-8 characters, while a genuine content loss differs by 60-260.
MATERIAL <- 25L

report <- list()
for (t in targets) {
    tmp <- tempfile(fileext = ".Rd")
    ok <- system2("git", c("-C", wd, "show", paste0("HEAD:", t$rd)),
                  stdout = tmp, stderr = FALSE) == 0L
    if (!ok) {
        report[[t$rd]] <- "could not read HEAD version"
        next
    }
    committed <- rd_args(tmp)
    unlink(tmp)
    if (!length(committed)) {
        report[[t$rd]] <- "no \\arguments in the committed Rd"
        next
    }
    if (is.null(t$fn)) {
        report[[t$rd]] <- "topic page - handled separately"
        next
    }
    src_path <- file.path(wd, t$src)
    lines <- readLines(src_path, warn = FALSE)
    def <- grep(paste0("(^|[^A-Za-z0-9._])", t$fn, "\\s*<-\\s*function"), lines)
    if (!length(def)) {
        report[[t$rd]] <- paste("function not found in", t$src)
        next
    }
    d <- def[[1L]]
    # The roxygen block is the run of comment lines immediately above the
    # definition. It must accept BOTH `#'` and plain `#` lines: roxygen puts
    # @importFrom on a plain `#` line, and stopping the walk there truncated
    # the block and hid every @param above it.
    start <- d
    while (start > 1L && grepl("^\\s*#", lines[start - 1L])) start <- start - 1L
    block <- lines[start:(d - 1L)]
    have <- block_params(block)

    missing <- setdiff(names(committed), names(have))
    thinner <- character()
    for (key in intersect(names(committed), names(have))) {
        src_txt <- have[[key]]$text
        rd_txt <- committed[[key]]$text
        if (nchar(rd_txt) - nchar(src_txt) >= MATERIAL) thinner <- c(thinner, key)
    }
    if (!length(missing) && !length(thinner)) {
        report[[t$rd]] <- "nothing missing"
        next
    }

    # 1) replace the @param entries whose text is thinner. A @param may span
    #    several `#'` lines; the whole entry has to be replaced, otherwise the
    #    continuation lines survive and the text ends up duplicated.
    out_block <- character()
    i <- 1L
    while (i <= length(block)) {
        m <- regmatches(block[i], regexec("^\\s*#'\\s*@param\\s+([^[:space:]]+)", block[i]))[[1]]
        if (!length(m)) {
            out_block <- c(out_block, block[i])
            i <- i + 1L
            next
        }
        key <- param_key(m[[2]])
        j <- i + 1L
        while (j <= length(block) && grepl("^\\s*#'", block[[j]]) &&
               !grepl("^\\s*#'\\s*@", block[[j]])) j <- j + 1L
        if (key %in% thinner) {
            out_block <- c(out_block, wrap_param(committed[[key]]$names, committed[[key]]$text))
        } else {
            out_block <- c(out_block, block[i:(j - 1L)])
        }
        i <- j
    }
    block <- out_block
    # 2) add the @param entries that are missing entirely
    new_lines <- character()
    for (key in missing) {
        new_lines <- c(new_lines, wrap_param(committed[[key]]$names, committed[[key]]$text))
    }
    report[[t$rd]] <- paste0(
        "add[", paste(missing, collapse = ", "), "] ",
        "reword[", paste(thinner, collapse = ", "), "]"
    )
    if (apply) {
        # insert after the last existing @param, else after @title
        param_idx <- grep("^\\s*#'\\s*@param\\s+", block)
        insert_at <- if (length(param_idx)) max(param_idx) + 1L else {
            ti <- grep("^\\s*#'\\s*@title", block)
            if (length(ti)) max(ti) + 1L else 1L
        }
        new_block <- append(block, new_lines, after = insert_at - 1L)
        out_lines <- c(lines[1:(start - 1L)], new_block, lines[d:length(lines)])
        writeLines(out_lines, src_path)
    }
}

cat(if (apply) "APPLIED\n" else "DRY RUN (re-run with --apply to write)\n\n")
for (rd in names(report)) cat(sprintf("  %-44s %s\n", rd, report[[rd]]))
