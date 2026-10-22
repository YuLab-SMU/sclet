#!/usr/bin/env Rscript
# audit-duplicate-definitions.R
#
# Finds every top-level name that is assigned more than once in the package:
#   (a) in two different R/ files, and
#   (b) twice inside one R/ file.
#
# Only names that end up in NAMESPACE (exported, or registered as an S3/S4
# method) are interesting, because a duplicate internal name is invisible to
# users -- but the report lists both so the "harmless" ones can be dismissed.
#
# Runtime resolution: R sources R/ files in DESCRIPTION's Collate order and
# later assignments overwrite earlier ones inside the package namespace.  The
# script therefore reports, for each duplicated name, which definition is last
# in Collate (the one that wins) and which is dead.
#
# Usage: Rscript .dev/tools/audit-duplicate-definitions.R [pkgdir]

pkg <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(pkg)) pkg <- "."
pkg <- normalizePath(pkg)
r_dir <- file.path(pkg, "R")

# ---- Collate order ---------------------------------------------------------
dcf <- read.dcf(file.path(pkg, "DESCRIPTION"))
collate_raw <- strsplit(trimws(sub("^Collate:\\s*", "", dcf[1, "Collate"])), "\\s+")[[1]]
collate <- gsub("['\"]", "", collate_raw)
collate <- collate[nzchar(collate)]
collate_rank <- stats::setNames(seq_along(collate), collate)

# ---- NAMESPACE exported names ----------------------------------------------
ns <- readLines(file.path(pkg, "NAMESPACE"))
ns <- ns[!grepl("^\\s*#", ns)]
exported <- unique(sub("^export\\((.*)\\)\\s*$", "\\1",
                        grep("^export\\(", ns, value = TRUE)))
exported <- gsub('"', "", exported)
s3 <- unique(sub("^S3method\\(([^,]+),.*$", "\\1", grep("^S3method\\(", ns, value = TRUE)))
s3 <- gsub('"', "", s3)
public <- unique(c(exported, s3))

# ---- scan top-level assignments -------------------------------------------
# A "definition" is a top-level `<name> <- ...` or `<name> = ...` assignment,
# or a top-level `function(...)` (rare, always assigned).  We ignore assignments
# inside function bodies -- those are locals.
top_level_names <- function(path) {
  exprs <- tryCatch(parse(path, keep.source = TRUE),
                    error = function(e) {
                      message("PARSE ERROR: ", path, ": ", conditionMessage(e))
                      NULL
                    })
  if (is.null(exprs)) return(character(0))
  out <- character(0)
  for (e in exprs) {
    if (is.call(e)) {
      op <- e[[1]]
      opn <- if (is.symbol(op)) as.character(op) else ""
      if (opn %in% c("<-", "=", "assign") && length(e) >= 3) {
        lhs <- e[[2]]
        if (is.symbol(lhs)) out <- c(out, as.character(lhs))
      } else if (opn %in% c("<<-") && length(e) >= 3) {
        lhs <- e[[2]]
        if (is.symbol(lhs)) out <- c(out, as.character(lhs))
      }
    }
  }
  out
}

files <- sort(list.files(r_dir, pattern = "[.]R$", full.names = TRUE))
defs <- list()   # name -> data.frame(file, line, col)
for (f in files) {
  exprs <- tryCatch(parse(f, keep.source = TRUE), error = function(e) NULL)
  if (is.null(exprs)) next
  src <- readLines(f, warn = FALSE)
  for (e in exprs) {
    if (!is.call(e)) next
    op <- e[[1]]
    opn <- if (is.symbol(op)) as.character(op) else ""
    if (!(opn %in% c("<-", "=") && length(e) >= 3)) next
    lhs <- e[[2]]
    if (!is.symbol(lhs)) next
    nm <- as.character(lhs)
    sr <- attr(e, "srcref")
    line <- if (!is.null(sr)) sr[[1]] else NA_integer_
    defs[[nm]] <- rbind(defs[[nm]],
                        data.frame(name = nm,
                                   file = basename(f),
                                   line = line,
                                   exported = nm %in% public,
                                   stringsAsFactors = FALSE))
  }
}

tab <- do.call(rbind, defs)
if (is.null(tab)) stop("no definitions found")

cat("=== ALL DUPLICATE TOP-LEVEL NAMES (any name) ===\n")
cnt <- table(tab$name)
dups <- names(cnt)[cnt > 1]
cat("total top-level names:", length(cnt), " duplicated:", length(dups), "\n\n")

tab$rank <- collate_rank[tab$file]
tab <- tab[order(tab$name, tab$rank), ]

report <- function(rows, label) {
  for (nm in unique(rows$name)) {
    r <- rows[rows$name == nm, ]
    if (nrow(r) < 2) next
    cat(sprintf("--- %s : %s ---\n", label, nm))
    for (i in seq_len(nrow(r))) {
      cat(sprintf("    %-28s %-40s collate#%-3s exported=%s\n",
                  r$file[i], paste0("line ", r$line[i]),
                  r$rank[i], r$exported[i]))
    }
    w <- r[which.max(r$rank), ]
    cat(sprintf("    >> WINS at runtime: %s:%s  (collate #%s)\n",
                w$file, w$line, w$rank))
    losers <- r[r$rank != w$rank, ]
    if (nrow(losers)) {
      for (i in seq_len(nrow(losers))) {
        cat(sprintf("    >> DEAD at runtime:  %s:%s\n",
                    losers$file[i], losers$line[i]))
      }
    }
    cat("\n")
  }
}

exp_rows <- tab[tab$exported, ]
report(exp_rows, "EXPORTED / PUBLIC")

cat("\n=== DUPLICATE INTERNAL NAMES (not exported) ===\n")
int_rows <- tab[!tab$exported, ]
report(int_rows, "INTERNAL")

# ---- names defined in more than one FILE -----------------------------------
cat("\n=== SUMMARY: exported names defined in >1 FILE ===\n")
e <- tab[tab$exported, ]
by <- split(e, e$name)
multi_file <- Filter(function(d) length(unique(d$file)) > 1, by)
cat("count:", length(multi_file), "\n")
for (nm in names(multi_file)) {
  d <- multi_file[[nm]]
  cat(sprintf("  %-32s files: %s\n", nm, paste(unique(d$file), collapse = ", ")))
}

cat("\n=== SUMMARY: exported names assigned TWICE IN ONE FILE ===\n")
multi_assign <- Filter(function(d) length(d$file) == 1 && nrow(d) > 1, by)
cat("count:", length(multi_assign), "\n")
for (nm in names(multi_assign)) {
  d <- multi_assign[[nm]]
  cat(sprintf("  %-32s %s lines %s\n", nm, d$file[1],
              paste(d$line, collapse = ", ")))
}

# ---- exported names with NO definition found (sanity) -----------------------
missing <- setdiff(public, unique(tab$name))
cat("\n=== public names with no plain top-level assignment found:", length(missing), "===\n")
if (length(missing)) cat(paste(missing, collapse = ", "), "\n")
