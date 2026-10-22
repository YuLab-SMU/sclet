# Run the \examples blocks of Rd files and report PASS / FAIL / SKIP.
#
#   Rscript .dev/tools/check-examples.R man/Foo.Rd man/Bar.Rd
#   Rscript .dev/tools/check-examples.R --group G2
#   Rscript .dev/tools/check-examples.R --all
#
# The package under test is the WORKTREE source, loaded once with
# devtools::load_all(). Each example is evaluated in its own fresh child
# environment, so an example cannot accidentally depend on a variable another
# example defined -- which is what R CMD check would do to it.
#
# \dontrun{} examples are reported SKIP: they are present but not executed.
# That is the intended shape for anything needing a Suggests package, network
# access or a Python environment.

args <- commandArgs(trailingOnly = TRUE)
wd <- getwd()
if (!length(args)) {
    cat("usage: Rscript .dev/tools/check-examples.R <rd>... | --group <G> | --all\n")
    quit(status = 2L)
}

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

is_dontrun <- function(rd) {
    grepl("\\dontrun", paste(readLines(rd, warn = FALSE), collapse = "\n"), fixed = TRUE)
}

# Strip \dontrun{...} blocks (and \donttest{...}) from example code, keeping
# the rest. A \dontrun body is never executed, so evaluating it would be
# wrong; but the setup lines that sit OUTSIDE it ARE executed by R CMD check,
# so they must be verified here too. Returns the code with the dontrun
# regions replaced by nothing, plus a flag telling whether anything is left.
strip_dontrun <- function(code) {
    out <- character()
    lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
    depth <- 0L
    for (ln in lines) {
        if (depth > 0L) {
            # count closing braces on this line to know when we are back out
            closes <- nchar(ln, type = "bytes") -
                nchar(gsub("}", "", ln, fixed = TRUE), type = "bytes")
            depth <- depth - closes
            if (depth < 0L) depth <- 0L
            next
        }
        if (grepl("\\\\(dontrun|donttest)\\s*\\{", ln)) {
            opens <- nchar(ln, type = "bytes") -
                nchar(gsub("{", "", ln, fixed = TRUE), type = "bytes")
            closes <- nchar(ln, type = "bytes") -
                nchar(gsub("}", "", ln, fixed = TRUE), type = "bytes")
            depth <- max(0L, opens - closes)
            next
        }
        out <- c(out, ln)
    }
    list(code = paste(out, collapse = "\n"), stripped = depth > 0L)
}

# Pull the @examples tag out of the roxygen block that documents `name` in
# `src_file`. This is the pre-roxygenise path: it lets an agent verify the
# example it just wrote without waiting for the (single, Lead-run)
# roxygenise() pass.
extract_examples_from_source <- function(src_file, name) {
    path <- file.path(wd, "R", src_file)
    if (!file.exists(path)) return(NULL)
    lines <- readLines(path, warn = FALSE)
    # locate the definition
    def <- grep(paste0("(^|[^A-Za-z0-9._])", name, "\\s*<-\\s*function"), lines)
    if (!length(def)) return(NULL)
    d <- def[[1L]]
    # walk back to the start of the contiguous roxygen block
    start <- d
    while (start > 1L && grepl("^#'", lines[start - 1L])) start <- start - 1L
    block <- lines[start:(d - 1L)]
    if (!length(block)) return(NULL)
    is_ex <- grep("^#'\\s*@examples", block)
    if (!length(is_ex)) return(NULL)
    i <- is_ex[[1L]]
    out <- character()
    # first line may carry inline code after the tag
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


if (identical(args[[1L]], "--all")) {
    files <- list.files(file.path(wd, "man"), pattern = "[.]Rd$", full.names = TRUE)
    mode <- "rd"
} else if (identical(args[[1L]], "--group")) {
    tsv <- file.path(wd, ".dev", "tools", "groups", paste0(args[[2L]], ".tsv"))
    if (!file.exists(tsv)) {
        cat("no such group file:", tsv, "\n")
        quit(status = 2L)
    }
    files <- file.path(wd, "man", sub("\t.*$", "", readLines(tsv)[-1]))
    mode <- "rd"
} else if (identical(args[[1L]], "--src-group")) {
    tsv <- file.path(wd, ".dev", "tools", "groups", paste0(args[[2L]], ".tsv"))
    if (!file.exists(tsv)) {
        cat("no such group file:", tsv, "\n")
        quit(status = 2L)
    }
    rows <- read.delim(tsv, stringsAsFactors = FALSE)
    files <- file.path(wd, "man", paste0(rows$name, ".Rd"))
    mode <- "src"
    src_rows <- rows
} else {
    files <- args
    mode <- "rd"
}

suppressWarnings(suppressMessages(devtools::load_all(wd, quiet = TRUE)))

status <- character(length(files))
notes <- character(length(files))
for (k in seq_along(files)) {
    rd <- files[[k]]
    nm <- basename(rd)
    if (identical(mode, "src")) {
        ex <- extract_examples_from_source(src_rows$src_file[[k]], src_rows$name[[k]])
        if (is.null(ex) || !nzchar(trimws(ex))) {
            status[k] <- "NOEXAMPLE"; next
        }
    } else {
        if (!file.exists(rd)) {
            status[k] <- "MISSING"; notes[k] <- "file not found"; next
        }
        ex <- extract_examples(rd)
        if (is.null(ex) || !nzchar(trimws(ex))) {
            status[k] <- "NOEXAMPLE"; next
        }
    }
    # A \dontrun body is never executed, but everything outside it is. So
    # strip the dontrun regions and evaluate what remains; only a page whose
    # ENTIRE example sits inside \dontrun has nothing left to run.
    st <- strip_dontrun(ex)
    runnable <- st$code
    if (!nzchar(trimws(runnable)) || !length(parse(text = runnable, keep.source = FALSE))) {
        status[k] <- "SKIP(dontrun)"; next
    }
    env <- new.env(parent = globalenv())
    res <- tryCatch({
        suppressWarnings(suppressMessages(eval(parse(text = runnable), envir = env)))
        "PASS"
    }, error = function(e) paste0("FAIL: ", conditionMessage(e)))
    status[k] <- res
    if (res != "PASS") {
        notes[k] <- paste(utils::head(utils::tail(trimws(runnable), 6), 3), collapse = " | ")
    }
}
res <- data.frame(rd = basename(files), status = status, note = notes, stringsAsFactors = FALSE)
print(res, row.names = FALSE)
cat("\nsummary:\n")
print(table(res$status))
if (any(res$status %in% c("FAIL", "NOEXAMPLE", "MISSING"))) quit(status = 1L)
