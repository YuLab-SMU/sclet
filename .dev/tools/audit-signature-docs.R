#!/usr/bin/env Rscript
# audit-signature-docs.R
#
# Cross-checks every exported function defined in R/ against its documentation.
#
# The "documented parameter" set is the Rd \arguments \item{} names (roxygen2
# turns the @param tags into exactly these).  The "signature" set is the
# function's real formals().  Functions that share one Rd topic (a getter and
# its setter, e.g. ActiveIdent / ActiveIdent<-) are compared against the union
# of their formals, because that is what the topic documents.
#
#   A. a formal argument with no \item       -> codoc "Undocumented arguments"
#   B. an \item for an argument that does not exist
#                                             -> codoc "Documented arguments not in \usage"
#   C. Rd \usage names an argument that is not a formal
#
# Usage: Rscript .dev/tools/audit-signature-docs.R [pkgdir]

pkg <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(pkg)) pkg <- "."
pkg <- normalizePath(pkg)

suppressMessages(devtools::load_all(pkg, quiet = TRUE))
ns <- asNamespace("sclet")

ns_lines <- readLines(file.path(pkg, "NAMESPACE"))
ns_lines <- ns_lines[!grepl("^\\s*#", ns_lines)]
exported <- unique(gsub('"', "", sub("^export\\((.*)\\)\\s*$", "\\1",
                                     grep("^export\\(", ns_lines, value = TRUE))))

# ---- index every Rd file by alias -------------------------------------------
rd_files <- list.files(file.path(pkg, "man"), pattern = "[.]Rd$", full.names = TRUE)

rd_block <- function(txt, what) {
  i <- grep(paste0("^\\\\", what, "\\s*\\{"), txt)
  if (!length(i)) return(character(0))
  j <- i[1]; depth <- 1L; k <- j + 1L; body <- character(0)
  while (k <= length(txt) && depth > 0) {
    depth <- depth + lengths(regmatches(txt[k], gregexpr("\\{", txt[k]))) -
      lengths(regmatches(txt[k], gregexpr("\\}", txt[k])))
    body <- c(body, txt[k]); k <- k + 1L
  }
  if (length(body)) body <- body[-length(body)]
  body
}

alias_index <- list()
rd_cache <- list()
for (f in rd_files) {
  txt <- readLines(f, warn = FALSE)
  al <- unlist(regmatches(txt, gregexpr("^\\\\alias\\{[^}]*\\}", txt)))
  al <- sub("^\\\\alias\\{", "", al); al <- sub("\\}$", "", al)
  for (a in al) if (is.null(alias_index[[a]])) alias_index[[a]] <- f
  args <- rd_block(txt, "arguments")
  items <- unlist(regmatches(args, gregexpr("\\\\item\\{[^}]*\\}", args)))
  items <- sub("^\\\\item\\{", "", items); items <- sub("\\}$", "", items)
  # an \item may name several arguments at once: \item{a, b}
  items <- unlist(strsplit(gsub("\\\\", "", items), ","))
  items <- trimws(items); items <- items[nzchar(items)]
  usage <- rd_block(txt, "usage")
  rd_cache[[basename(f)]] <- list(items = unique(items), usage = usage, alias = al)
}

# ---- usage: argument names of the call in \usage -----------------------------
usage_calls <- function(usage_lines) {
  if (!length(usage_lines)) return(list())
  txt <- paste(usage_lines, collapse = "\n")
  out <- list()
  pos <- 1L
  while (TRUE) {
    m <- regexpr("(^|[^A-Za-z0-9._])([A-Za-z._][A-Za-z0-9._]*)\\s*\\(", substr(txt, pos, nchar(txt)), perl = TRUE)
    if (m < 0) break
    st <- pos + attr(m, "capture.start")[2] - 1L
    en <- pos + attr(m, "capture.start")[2] - 1L
    # find the matching close paren
    depth <- 0L; k <- en; chars <- strsplit(substr(txt, en, nchar(txt)), "")[[1]]
    inner <- character(0)
    for (ch in chars) {
      if (ch == "(") depth <- depth + 1L
      if (ch == ")") { depth <- depth - 1L; if (depth == 0L) break }
      inner <- c(inner, ch)
    }
    inner_txt <- paste(inner, collapse = "")
    parts <- strsplit(inner_txt, ",", fixed = TRUE)[[1]]
    parts <- trimws(parts); parts <- parts[nzchar(parts)]
    nm <- trimws(sub("=.*$", "", parts))
    out[[length(out) + 1L]] <- list(name = substr(txt, st, en - 1L), args = nm)
    pos <- en + length(inner) + 1L
    if (pos > nchar(txt)) break
  }
  out
}

# ---- gather ------------------------------------------------------------------
info <- list()
for (nm in exported) {
  if (!exists(nm, envir = ns, inherits = FALSE)) next
  obj <- get(nm, envir = ns, inherits = FALSE)
  if (!is.function(obj)) next
  fa <- names(formals(obj))
  if (is.null(fa)) next
  rdf <- alias_index[[nm]]
  info[[nm]] <- list(name = nm, formals = fa,
                     rd = if (is.null(rdf)) NA_character_ else basename(rdf))
}

# group by Rd topic
by_rd <- list()
for (nm in names(info)) {
  r <- info[[nm]]$rd
  if (is.na(r)) next
  by_rd[[r]] <- c(by_rd[[r]], nm)
}

rows <- list()
for (r in names(by_rd)) {
  nms <- by_rd[[r]]
  union_formals <- unique(unlist(lapply(nms, function(n) info[[n]]$formals)))
  items <- rd_cache[[r]]$items
  u <- usage_calls(rd_cache[[r]]$usage)
  # for a topic with several functions, keep the union of the documented calls
  uargs <- unique(unlist(lapply(u, function(x) x$args)))
  rows[[length(rows) + 1L]] <- data.frame(
    rd = r,
    functions = paste(nms, collapse = " / "),
    union_formals = paste(union_formals, collapse = ","),
    rd_items = paste(items, collapse = ","),
    usage_args = paste(uargs, collapse = ","),
    undocumented = paste(setdiff(setdiff(union_formals, "..."), items), collapse = ","),
    extra_doc = paste(setdiff(items, c(union_formals, "...")), collapse = ","),
    usage_mismatch = paste(setdiff(setdiff(union_formals, "..."), uargs), collapse = ","),
    usage_extra = paste(setdiff(uargs, union_formals), collapse = ","),
    stringsAsFactors = FALSE
  )
}
res <- do.call(rbind, rows)
cat("Rd topics covering exported functions:", nrow(res), "\n\n")

show <- function(d, label, col) {
  cat(sprintf("\n=== %s ===\ncount: %d\n", label, nrow(d)))
  if (!nrow(d)) return(invisible())
  for (i in seq_len(nrow(d))) {
    cat(sprintf("  %-36s %-30s %s\n", d$functions[i], d$rd[i], d[[col]][i]))
  }
}

show(res[res$undocumented != "", , drop = FALSE],
     "A. formal argument with NO \\item in the Rd \\arguments section", "undocumented")
show(res[res$extra_doc != "", , drop = FALSE],
     "B. Rd \\arguments documents an argument that does not exist", "extra_doc")
show(res[res$usage_mismatch != "", , drop = FALSE],
     "C. Rd \\usage omits a formal argument", "usage_mismatch")
show(res[res$usage_extra != "", , drop = FALSE],
     "D. Rd \\usage names an argument that is not a formal", "usage_extra")

saveRDS(res, file.path(pkg, ".dev", "tools", ".audit-signature-docs.rds"))

# ---- also dump every exported signature, for the drift analysis --------------
sig <- do.call(rbind, lapply(names(info), function(nm) {
  data.frame(name = nm,
             file = NA_character_, line = NA_integer_,
             signature = paste(info[[nm]]$formals, collapse = ", "),
             stringsAsFactors = FALSE)
}))
# attach source file:line
srcref_map <- list()
for (f in list.files(file.path(pkg, "R"), pattern = "[.]R$", full.names = TRUE)) {
  exprs <- tryCatch(parse(f, keep.source = TRUE), error = function(e) NULL)
  if (is.null(exprs)) next
  srefs <- attr(exprs, "srcref")
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (!is.call(e)) next
    op <- e[[1]]
    if (!(is.symbol(op) && as.character(op) %in% c("<-", "=") && length(e) >= 3)) next
    lhs <- e[[2]]
    if (!is.symbol(lhs)) next
    sr <- if (is.null(srefs)) NULL else srefs[[i]]
    srcref_map[[as.character(lhs)]] <- c(basename(f), if (is.null(sr)) NA_integer_ else sr[[1]])
  }
}
for (i in seq_len(nrow(sig))) {
  m <- srcref_map[[sig$name[i]]]
  if (!is.null(m)) { sig$file[i] <- m[1]; sig$line[i] <- as.integer(m[2]) }
}
saveRDS(sig, file.path(pkg, ".dev", "tools", ".audit-signatures.rds"))
cat("\nwrote .dev/tools/.audit-signature-docs.rds and .audit-signatures.rds\n")
