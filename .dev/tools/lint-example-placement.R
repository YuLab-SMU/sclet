# Verify that every @examples block sits in the file that actually DEFINES
# the function it documents.
#
# The group partition was built from a name->file lookup, and that lookup was
# wrong for a while (an Rd file name is not an object name: subset-cell.Rd
# documents subset_cell; a @rdname topic page documents several functions).
# This check does not trust the partition -- it re-derives the answer from the
# source and compares.
#
#   Rscript .dev/tools/lint-example-placement.R

wd <- getwd()

# which file defines each top-level assignment
defs <- list()
for (f in list.files(file.path(wd, "R"), pattern = "[.]R$", full.names = TRUE)) {
    exprs <- tryCatch(parse(f, keep.source = FALSE), error = function(e) NULL)
    if (is.null(exprs)) next
    for (e in exprs) {
        if (is.call(e) && length(e) >= 3L &&
            identical(as.character(e[[1]]), "<-") && is.name(e[[2]])) {
            defs[[as.character(e[[2]])]] <- basename(f)
        }
    }
}

problems <- list()
total <- 0L
for (f in list.files(file.path(wd, "R"), pattern = "[.]R$", full.names = TRUE)) {
    lines <- readLines(f, warn = FALSE)
    base <- basename(f)
    i <- 1L
    while (i <= length(lines)) {
        # find a roxygen block start
        if (!grepl("^\\s*#'", lines[i])) {
            i <- i + 1L
            next
        }
        start <- i
        while (i <= length(lines) && grepl("^\\s*#", lines[i])) i <- i + 1L
        block <- lines[start:(i - 1L)]
        # what follows the block?
        j <- i
        while (j <= length(lines) && !nzchar(trimws(lines[j]))) j <- j + 1L
        if (j > length(lines)) next
        target <- sub("^\\s*([^\\s<-]+)\\s*<-.*$", "\\1", lines[j])
        if (identical(target, lines[j])) {
            i <- j
            next
        }
        if (target %in% defs && !identical(defs[[target]], base)) {
            problems[[paste0(base, ":", target)]] <- defs[[target]]
        }
        i <- j
    }
    total <- total + 1L
}

# now the @examples-specific check: the block that carries @examples must
# document a function defined in the same file
misplaced <- list()
for (f in list.files(file.path(wd, "R"), pattern = "[.]R$", full.names = TRUE)) {
    lines <- readLines(f, warn = FALSE)
    base <- basename(f)
    i <- 1L
    while (i <= length(lines)) {
        if (!grepl("^\\s*#'", lines[i])) {
            i <- i + 1L
            next
        }
        start <- i
        while (i <= length(lines) && grepl("^\\s*#", lines[i])) i <- i + 1L
        block <- lines[start:(i - 1L)]
        if (!any(grepl("@examples", block))) next
        j <- i
        while (j <= length(lines) && !nzchar(trimws(lines[j]))) j <- j + 1L
        if (j > length(lines)) next
        target <- sub("^\\s*([^\\s<-]+)\\s*<-.*$", "\\1", lines[j])
        if (identical(target, lines[j])) next
        if (target %in% defs && !identical(defs[[target]], base)) {
            misplaced[[paste0(base, " -> ", target)]] <- defs[[target]]
        }
        i <- j
    }
}

cat("roxygen blocks whose function is defined elsewhere (all blocks):\n")
if (length(problems)) {
    for (k in names(problems)) cat(sprintf("  %-52s defined in %s\n", k, problems[[k]]))
} else {
    cat("  none\n")
}
cat("\n@examples blocks placed in the wrong file:\n")
if (length(misplaced)) {
    for (k in names(misplaced)) cat(sprintf("  %-52s defined in %s\n", k, misplaced[[k]]))
} else {
    cat("  none\n")
}
if (length(misplaced)) quit(status = 1L)
