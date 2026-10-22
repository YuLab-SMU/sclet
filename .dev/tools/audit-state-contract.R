#!/usr/bin/env Rscript
# audit-state-contract.R
#
# Measures, per analysis module, how many exported functions in that module
# register an analysis-state record (sclet_set_analysis_state()) and how many
# write their result onto the object without one.
#
# Also records, per exported function:
#   - which check_*_readiness() gate it calls (if any)
#   - whether it stops on error (stop()/cli::cli_abort())
#   - whether it swallows errors (tryCatch with a handler that returns
#     NULL / a default)
#   - whether it calls sclet_log_command()
#
# Usage: Rscript .dev/tools/audit-state-contract.R [pkgdir] > report.txt

pkg <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(pkg)) pkg <- "."
pkg <- normalizePath(pkg)

ns <- asNamespace("sclet")
src <- file.path(pkg, "R")

# ---- collect every top-level function with its file, line and body ----------
recs <- list()
for (f in list.files(src, pattern = "[.]R$", full.names = TRUE)) {
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
    nm <- as.character(lhs)
    sr <- if (is.null(srefs)) NULL else srefs[[i]]
    body_txt <- paste(deparse(e[[3]]), collapse = "\n")
    recs[[nm]] <- list(name = nm, file = basename(f),
                       line = if (is.null(sr)) NA_integer_ else sr[[1]],
                       body = body_txt)
  }
}

ns_lines <- readLines(file.path(pkg, "NAMESPACE"))
ns_lines <- ns_lines[!grepl("^\\s*#", ns_lines)]
exported <- unique(gsub('"', "", sub("^export\\((.*)\\)\\s*$", "\\1",
                                     grep("^export\\(", ns_lines, value = TRUE))))

gate_fns <- sort(grep("^check_.*readiness$", exported, value = TRUE))

# ---- per-function facts -----------------------------------------------------
rows <- list()
for (nm in exported) {
  r <- recs[[nm]]
  if (is.null(r)) next
  b <- r$body
  gates <- gate_fns[vapply(gate_fns, function(g) grepl(paste0("\\b", g, "\\s*\\("), b), logical(1))]
  rows[[length(rows) + 1L]] <- data.frame(
    name = nm,
    file = r$file,
    line = r$line,
    core = !grepl("^ai-", r$file),
    sets_state = grepl("sclet_set_analysis_state", b),
    sets_record = grepl("sclet_set_state_record", b),
    sets_any = grepl("\\bsclet_set_[A-Za-z0-9_]+\\(", b),
    logs_command = grepl("\\bsclet_log_command\\(", b),
    stops = grepl("\\bstop\\(|cli::cli_abort\\(|abort\\(", b),
    trycatch = grepl("\\btryCatch\\(", b),
    gates = paste(gates, collapse = ","),
    stringsAsFactors = FALSE
  )
}
res <- do.call(rbind, rows)

core <- res[res$core & grepl("^(Run|Find|Call|Classify|Aggregate|Compare|Compute|Score|Detect|Batch|Normalize|Scale|Subset|subset|Merge|Write|Read)", res$name), ]
core <- core[!is.na(core$gates) | TRUE, ]

cat("################ CORE (non-ai) exported functions analysed ################\n\n")

# ---- module-level table -----------------------------------------------------
core$module <- sub("\\.R$", "", core$file)
mods <- sort(unique(core$module))
cat(sprintf("%-28s %5s %5s %5s %5s %5s %5s\n",
            "module", "n", "state", "record", "noState", "gate", "noGate"))
for (m in mods) {
  d <- core[core$module == m, ]
  n <- nrow(d)
  st <- sum(d$sets_state)
  sr <- sum(d$sets_record)
  gt <- sum(d$gates != "")
  cat(sprintf("%-28s %5d %5d %5d %5d %5d %5d\n", m, n, st, sr, n - sum(d$sets_state | d$sets_record), gt, n - gt))
}
cat("\nTOTALS  n=", nrow(core),
    " setsAnalysisState=", sum(core$sets_state),
    " setsStateRecord=", sum(core$sets_record),
    " either=", sum(core$sets_state | core$sets_record),
    " neither=", sum(!core$sets_state & !core$sets_record),
    " callsGate=", sum(core$gates != ""),
    " noGate=", sum(core$gates == ""), "\n", sep = "")

cat("\n################ CORE exported functions WITHOUT a state record ################\n")
d <- core[!core$sets_state & !core$sets_record, ]
d <- d[order(d$file, d$line), ]
for (i in seq_len(nrow(d))) {
  cat(sprintf("  %-34s %-26s %s:%s   gates=[%s] stops=%s writes=%s\n",
              d$name[i], d$file[i], d$file[i], d$line[i], d$gates[i],
              d$stops[i], d$sets_any[i]))
}

cat("\n################ CORE exported functions WITH a state record ################\n")
d <- core[core$sets_state, ]
d <- d[order(d$file, d$line), ]
for (i in seq_len(nrow(d))) {
  cat(sprintf("  %-34s %-26s %s:%s   gates=[%s]\n",
              d$name[i], d$file[i], d$file[i], d$line[i], d$gates[i]))
}

cat("\n################ GATES: which core functions call which gate ################\n")
for (g in gate_fns) {
  callers <- core[grepl(paste0("\\b", g, "\\b"), core$gates), ]
  cat(sprintf("\n%s -- called by %d core function(s)\n", g, nrow(callers)))
  for (i in seq_len(nrow(callers))) {
    cat(sprintf("    %-34s %s:%s\n", callers$name[i], callers$file[i], callers$line[i]))
  }
}

cat("\n################ gate functions: are they exported? ################\n")
for (g in gate_fns) {
  r <- recs[[g]]
  if (is.null(r)) { cat(g, ": NOT FOUND in R/\n"); next }
  cat(sprintf("  %-38s %s:%s\n", g, r$file, r$line))
}

saveRDS(res, file.path(pkg, ".dev", "tools", ".audit-state-contract.rds"))
cat("\nwrote .dev/tools/.audit-state-contract.rds\n")
