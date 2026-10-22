#!/usr/bin/env Rscript
# Record what a bookdown render actually EXECUTED, and which sclet built it.
#
# Why this exists: the book executes ~194 R chunks against real data, and it is
# the only place many analysis entry points are exercised at all (69 of the 283
# exported functions appear in the book and nowhere in tests). The book is
# rendered locally on purpose -- the Python-backed steps make a CI run far too
# heavy -- which means nothing in CI would notice if a chapter stopped
# executing. Without a record, "it ran last time" is a memory, not a fact.
#
# So after a render this writes `_book_full/render-record.md`: the sclet commit,
# the R session, and per chapter the number of chunks, how many were executed,
# and how many were skipped. A chapter whose executed count drops to zero is
# visible by diffing two records, which is the only reliable way to see a
# chapter that was quietly disabled.
#
# Usage (after `make gitbook`):
#   Rscript record-render.R [book_dir] [output_file]

args <- commandArgs(trailingOnly = TRUE)
book_dir <- if (length(args) >= 1L) args[[1L]] else "."
out_file <- if (length(args) >= 2L) args[[2L]] else "_book_full/render-record.md"

book_dir <- normalizePath(book_dir, mustWork = TRUE)
rmd_files <- list.files(book_dir, pattern = "[.]Rmd$", full.names = TRUE)
if (!length(rmd_files)) {
    stop("no .Rmd files in ", book_dir, call. = FALSE)
}

# The chunk header line: ```{r ...}  or  ```{r ... eval=FALSE}
chunk_re <- "^```[:space:]*\\{r\\b([^}]*)\\}"
chunk_end_re <- "^```[:space:]*$"

# Chunks fall into three groups, and the difference between the last two matters:
#
#   executed  -- no `eval=` option, or `eval=TRUE`. Always runs.
#   default   -- `eval=` is a VARIABLE. Whether it runs depends on the variable,
#                and the book defines those variables with defaults:
#                  eval_chunks  <- SCLET_BOOKDOWN_EVAL, default TRUE
#                  eval_ai      <- ... EVAL_AI=1 and aisdk, default FALSE
#                  svp_ready / pipeline_ready / *_cache_ready / file.exists /
#                  exists / requireNamespace -- data-dependent, false on a
#                  machine without that data.
#                Counting these as "executed" would overstate the record (the AI
#                chapters would look run when they never ran); counting them as
#                "off" would understate it (eval_chunks is TRUE by default, so
#                65 chunks DO run). So they are recorded as `conditional`, and
#                a run should be judged by the numbers together, not by one of
#                them.
#   off       -- `eval=FALSE`. Never runs.
ai_switch_on <- identical(
    tolower(Sys.getenv("SCLET_BOOKDOWN_EVAL_AI", "false")),
    "true"
) || identical(Sys.getenv("SCLET_BOOKDOWN_EVAL_AI"), "1") ||
    identical(tolower(Sys.getenv("SCLET_BOOKDOWN_EVAL_AI")), "yes")

scan_file <- function(path) {
    lines <- readLines(path, warn = FALSE)
    header <- grep(chunk_re, lines)
    if (!length(header)) {
        return(list(chunks = 0L, executed = 0L, off = 0L, conditional = 0L,
                    ai_gated = 0L, ai_ran = 0L))
    }
    opts <- sub(chunk_re, "\\1", lines[header])
    is_off <- grepl("eval\\s*=\\s*FALSE", opts)
    is_var <- grepl("eval\\s*=\\s*[A-Za-z_][A-Za-z0-9_.]*\\s*$", opts) & !is_off
    is_ai <- grepl("eval\\s*=\\s*eval_ai", opts)
    # `eval_ai` is the ONE switch whose value this script can know for certain:
    # it is `EVAL_AI=1 AND aisdk installed AND eval_chunks`. `eval_chunks`
    # (SCLET_BOOKDOWN_EVAL) defaults to TRUE, so it does not subtract. So when
    # the switch was on, those chunks RAN, and calling them "conditional" would
    # make a run with the AI chapters enabled look identical to one without --
    # which is exactly the distinction the record exists to capture.
    ai_ran <- is_ai & ai_switch_on
    list(
        chunks = length(header),
        executed = sum((!is_off & !is_var) | ai_ran),
        off = sum(is_off),
        conditional = sum(is_var & !ai_ran),
        ai_gated = sum(is_ai),
        ai_ran = sum(ai_ran)
    )
}

sclet_commit <- function(dir) {
    out <- tryCatch(
        suppressWarnings(system2(
            "git", c("-C", dir, "rev-parse", "--short", "HEAD"),
            stdout = TRUE, stderr = FALSE
        )),
        error = function(e) character()
    )
    if (length(out) == 1L && nzchar(out)) out else NA_character_
}

sclet_version <- tryCatch(
    as.character(utils::packageVersion("sclet")),
    error = function(e) NA_character_
)
installed_from <- tryCatch(
    dirname(dirname(system.file(package = "sclet"))),
    error = function(e) NA_character_
)

rows <- lapply(sort(basename(rmd_files)), function(f) {
    s <- scan_file(file.path(book_dir, f))
    data.frame(file = f, chunks = s$chunks, executed = s$executed,
               conditional = s$conditional, ai_gated = s$ai_gated,
               ai_ran = s$ai_ran, off = s$off, stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)

# An .Rmd with chunks but zero executed is the failure shape this record
# exists to surface: it renders, so the build looks green, and every result on
# the page is stale prose.
dead <- tab$file[tab$chunks > 0L & tab$executed == 0L & tab$conditional == 0L]

lines <- c(
    "# Book render record",
    "",
    sprintf("- Rendered: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("- sclet version: %s", sclet_version),
    sprintf("- sclet library: %s", installed_from),
    sprintf("- sclet commit: %s", sclet_commit(book_dir)),
    sprintf("- R: %s", R.version.string),
    sprintf("- Platform: %s", R.version$platform),
    # Whether the AI chunks actually ran on this render. Without this line a
    # record cannot be read: the same source produces two different runs
    # depending on this switch, and the difference is exactly the part that
    # cannot be checked in CI.
    sprintf("- SCLET_BOOKDOWN_EVAL: %s",
            if (nzchar(Sys.getenv("SCLET_BOOKDOWN_EVAL")))
                Sys.getenv("SCLET_BOOKDOWN_EVAL") else "(unset, defaults to true)"),
    sprintf("- SCLET_BOOKDOWN_EVAL_AI: %s",
            if (nzchar(Sys.getenv("SCLET_BOOKDOWN_EVAL_AI")))
                Sys.getenv("SCLET_BOOKDOWN_EVAL_AI") else "(unset, AI chunks skipped)"),
    sprintf("- aisdk available: %s",
            requireNamespace("aisdk", quietly = TRUE)),
    "",
    "Counts are per chapter, over three groups:",
    "",
    "- **executed** — no `eval=` option, or `eval=TRUE`. Always runs.",
    "- **conditional** — `eval=` is a variable. `eval_chunks` (the",
    "  `SCLET_BOOKDOWN_EVAL` switch) defaults to TRUE, so most of these DO run;",
    "  `eval_ai` defaults to FALSE; `svp_ready` / `pipeline_ready` /",
    "  `file.exists` / `requireNamespace` depend on data this machine may not",
    "  have.",
    "- **off** — `eval=FALSE`. Never runs.",
    "",
    "**A record is only comparable to another record taken on a machine with",
    "the same data.** The point of the record is the diff: a chapter whose",
    "numbers drop means something stopped running.",
    "",
    "| chapter | chunks | executed | conditional | AI-gated (ran?) | off |",
    "|---|---|---|---|---|---|",
    sprintf("| %s | %d | %d | %d | %d (%s) | %d |",
            tab$file, tab$chunks, tab$executed, tab$conditional,
            tab$ai_gated,
            ifelse(tab$ai_gated == 0L, "-",
                   ifelse(tab$ai_ran == tab$ai_gated, "ran",
                          ifelse(tab$ai_ran == 0L, "skipped", "partial"))),
            tab$off),
    "",
    sprintf("**Totals**: %d chunks, %d executed, %d conditional, %d off; of the %d AI-gated chunks, %d ran.",
            sum(tab$chunks), sum(tab$executed), sum(tab$conditional),
            sum(tab$off), sum(tab$ai_gated), sum(tab$ai_ran)),
    ""
)

if (length(dead)) {
    lines <- c(lines, "## Chapters with chunks but nothing executed", "",
               paste0("- ", dead), "")
}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, out_file)
cat(paste(lines, collapse = "\n"), "\n")
if (length(dead)) {
    cat("\nWARNING: ", length(dead), " chapter(s) executed no code: ",
        paste(dead, collapse = ", "), "\n", sep = "")
}
