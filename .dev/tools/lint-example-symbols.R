# Lint the @examples blocks for symbols that are never defined.
#
# A \dontrun{} body is never executed, so a broken call inside it is invisible
# to the checker -- and to R CMD check. But a reader who copies the example
# gets an error, and an example that cannot work teaches the wrong thing.
#
# This walks each example, collects the names it binds (assignments, function
# arguments, function definitions) and reports the names it reads that it never
# binds and that are not base R / package exports.
#
#   Rscript .dev/tools/lint-example-symbols.R --src-group G3
#   Rscript .dev/tools/lint-example-symbols.R --all

args <- commandArgs(trailingOnly = TRUE)
wd <- getwd()
if (!length(args)) {
    cat("usage: Rscript .dev/tools/lint-example-symbols.R <rd>... | --group <G> | --src-group <G> | --all\n")
    quit(status = 2L)
}

suppressWarnings(suppressMessages(devtools::load_all(wd, quiet = TRUE)))

base_names <- c(
    ls(baseenv()), ls(globalenv()),
    # the package's own exports
    getNamespaceExports("sclet")
)
# common free-text helpers that are fine in an example
allowed <- c("x", "i", "n")

extract_examples <- function(rd) {
    lines <- readLines(rd, warn = FALSE)
    start <- which(startsWith(lines, "\\examples{"))
    if (!length(start)) return(NULL)
    i <- start[[1L]]
    # The block ends at the `}` that closes \examples{ -- NOT at the first
    # `}` at the start of a line: a \dontrun{} inside the example closes
    # with a `}` in column 1 as well, and stopping there truncates the
    # example (which then fails to parse, and silently hides everything
    # after the \dontrun block from being checked).
    depth <- 1L
    end <- i
    while (end < length(lines)) {
        end <- end + 1L
        opens <- nchar(lines[end], type = "bytes") -
            nchar(gsub("{", "", lines[end], fixed = TRUE), type = "bytes")
        closes <- nchar(lines[end], type = "bytes") -
            nchar(gsub("}", "", lines[end], fixed = TRUE), type = "bytes")
        depth <- depth + opens - closes
        if (depth <= 0L) break
    }
    if (end > length(lines)) return(NULL)
    body <- sub("^    ", "", lines[(i + 1L):(end - 1L)])
    txt <- paste(body, collapse = "\n")
    # Roxygen escapes a literal % as \% in the Rd. Reading the Rd as plain
    # text and parsing it then fails on "invalid \% escape" -- so undo the
    # escaping before handing the code to parse().
    txt <- gsub("\\\\%", "%", txt)
    txt
}

extract_examples_from_source <- function(src_file, name) {
    path <- file.path(wd, "R", src_file)
    if (!file.exists(path)) return(NULL)
    lines <- readLines(path, warn = FALSE)
    def <- grep(paste0("(^|[^A-Za-z0-9._])", name, "\\s*<-\\s*function"), lines)
    if (!length(def)) return(NULL)
    d <- def[[1L]]
    start <- d
    while (start > 1L && grepl("^#'", lines[start - 1L])) start <- start - 1L
    block <- lines[start:(d - 1L)]
    if (!length(block)) return(NULL)
    is_ex <- grep("^#'\\s*@examples", block)
    if (!length(is_ex)) return(NULL)
    i <- is_ex[[1L]]
    out <- character()
    first <- sub("^#'\\s*@examples\\s*", "", block[i])
    if (nzchar(first)) out <- c(out, first)
    j <- i + 1L
    while (j <= length(block) && grepl("^#'", block[j])) {
        ln <- sub("^#'\\s?", "", block[j])
        if (grepl("^@", ln)) break
        out <- c(out, ln)
        j <- j + 1L
    }
    paste(out, collapse = "\n")
}

# Every name bound anywhere in the example (including inside dontrun)
#
# NOTE: "\dontrun{" is not valid R, so parse() FAILS on any example that
# contains one -- and a failed parse would silently make both walkers return
# empty, i.e. report "no undefined symbols" for every dontrun example. That
# blind spot hid five real defects. Normalise the marker to `if (FALSE) {`
# before parsing: the body still parses (so its names are seen) but it is
# never executed.
normalise_for_parse <- function(code) {
    code <- gsub("\\\\dontrun\\s*\\{", "if (FALSE) {", code)
    code <- gsub("\\\\donttest\\s*\\{", "if (FALSE) {", code)
    code
}

bound_names <- function(code) {
    exprs <- tryCatch(parse(text = normalise_for_parse(code), keep.source = FALSE),
                       error = function(e) e)
    if (inherits(exprs, "error")) {
        cat("  (parse failed after normalisation -- names not checked): ",
            conditionMessage(exprs), "\n")
        return(character())
    }
    nms <- character()
    walk <- function(e) {
        if (is.symbol(e)) {
            nms <<- c(nms, as.character(e))
        } else if (is.call(e)) {
            if (identical(as.character(e[[1]]), "<-") && is.name(e[[2]])) {
                nms <<- c(nms, as.character(e[[2]]))
            }
            if (identical(as.character(e[[1]]), "=") && is.name(e[[2]])) {
                nms <<- c(nms, as.character(e[[2]]))
            }
            if (identical(as.character(e[[1]]), "for") && is.name(e[[2]])) {
                nms <<- c(nms, as.character(e[[2]]))
            }
            if (identical(as.character(e[[1]]), "function")) {
                f <- e[[2]]
                if (is.name(f)) nms <<- c(nms, as.character(f))
                if (is.call(f) && identical(as.character(f[[1]]), "function")) {
                    for (a in as.list(f[[2]])) if (is.name(a)) nms <<- c(nms, as.character(a))
                }
            }
            for (a in as.list(e)) if (!missing(a)) walk(a)
        }
    }
    for (e in exprs) walk(e)
    unique(nms)
}

# every name READ as a value (not in call position, not a function name)
read_names <- function(code) {
    exprs <- tryCatch(parse(text = normalise_for_parse(code), keep.source = FALSE),
                       error = function(e) NULL)
    if (is.null(exprs)) return(character())
    nms <- character()
    walk <- function(e, in_call_head = FALSE) {
        if (is.symbol(e)) {
            if (!in_call_head) nms <<- c(nms, as.character(e))
            return()
        }
        if (is.call(e)) {
            head <- e[[1L]]
            walk(head, TRUE)
            for (a in as.list(e)[-1]) if (!missing(a)) walk(a, FALSE)
            return()
        }
        if (is.pairlist(e) || is.expression(e)) {
            for (a in as.list(e)) if (!missing(a)) walk(a, FALSE)
        }
    }
    for (e in exprs) walk(e, FALSE)
    unique(nms)
}

if (identical(args[[1L]], "--all")) {
    files <- list.files(file.path(wd, "man"), pattern = "[.]Rd$", full.names = TRUE)
    mode <- "rd"
} else if (identical(args[[1L]], "--group") || identical(args[[1L]], "--src-group")) {
    tsv <- file.path(wd, ".dev", "tools", "groups", paste0(args[[2L]], ".tsv"))
    rows <- read.delim(tsv, stringsAsFactors = FALSE)
    files <- file.path(wd, "man", paste0(rows$name, ".Rd"))
    if (identical(args[[1L]], "--src-group")) {
        mode <- "src"; src_rows <- rows
    } else {
        mode <- "rd"
    }
} else {
    files <- args
    mode <- "rd"
}

bad <- list()
for (k in seq_along(files)) {
    rd <- files[[k]]
    nm <- basename(rd)
    if (identical(mode, "src")) {
        ex <- extract_examples_from_source(src_rows$src_file[[k]], src_rows$name[[k]])
    } else {
        ex <- extract_examples(rd)
    }
    if (is.null(ex) || !nzchar(trimws(ex))) next
    b <- bound_names(ex)
    r <- read_names(ex)
    unknown <- setdiff(r, c(b, base_names, allowed))
    if (length(unknown)) {
        bad[[nm]] <- unknown
    }
}
if (length(bad)) {
    cat("examples reading names they never bind:\n\n")
    for (nm in names(bad)) cat(sprintf("  %-42s %s\n", nm, paste(bad[[nm]], collapse = ", ")))
    cat("\ntotal:", length(bad), "\n")
    quit(status = 1L)
}
cat("no undefined symbols found in", length(files), "examples\n")
