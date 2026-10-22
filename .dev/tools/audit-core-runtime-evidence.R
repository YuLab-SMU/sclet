#!/usr/bin/env Rscript
# audit-core-runtime-evidence.R
#
# Re-runnable evidence for the core-layer audit findings.  Every claim in
# .dev/audit-core.md that needs a runtime check is produced here.
#
#   E1  exported accessors that return a DIFFERENT record when asked for an
#       id that does not exist
#   E2  the "argument is of length zero" error path in the CCI group resolver
#   E3  the readiness-gate call graph (which gate is called from where)
#   E4  the package each readiness gate checks vs the package the analysis
#       actually requires
#   E5  RunCCI()'s keep_object argument reaching RunNicheNet(), which has no
#       such argument and no `...`
#
# Usage: Rscript .dev/tools/audit-core-runtime-evidence.R [pkgdir]

pkg <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(pkg)) pkg <- "."
pkg <- normalizePath(pkg)
suppressMessages(devtools::load_all(pkg, quiet = TRUE))

cat("=========== E1  accessor returns the WRONG record for a bogus id ===========\n")
types <- c("trajectory", "milo", "supercell", "cellchat", "cellphonedb", "nichenet",
           "velocity", "scenic", "geneset_scoring", "cellrank", "colocalization",
           "spatial_niche", "rare_cells", "augur", "celloracle", "batch")
sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
recs <- list(); acts <- list()
for (t in types) {
  recs[[t]] <- list(real_record = list(id = paste0(t, "-real"), type = t, method = "real"))
  acts[[t]] <- list(id = paste0(t, "-legacy"), marker = "LEGACY_SLOT")
}
S4Vectors::metadata(sce)$sclet <- list(
  states = list(active = lapply(types, function(t) paste0(t, "-real")), records = recs),
  analyses = acts)

accs <- list(
  get_trajectory      = function(o) get_trajectory(o, id = "BOGUS"),
  get_milo            = function(o) get_milo(o, id = "BOGUS"),
  get_supercell       = function(o) get_supercell(o, id = "BOGUS"),
  get_cellchat        = function(o) get_cellchat(o, id = "BOGUS"),
  get_cellphonedb     = function(o) get_cellphonedb(o, id = "BOGUS"),
  get_nichenet        = function(o) get_nichenet(o, id = "BOGUS"),
  get_velocity        = function(o) get_velocity(o, id = "BOGUS"),
  get_scenic          = function(o) get_scenic(o, id = "BOGUS"),
  get_geneset_scoring = function(o) get_geneset_scoring(o, id = "BOGUS"),
  get_cellrank        = function(o) get_cellrank(o, id = "BOGUS"),
  get_spatial         = function(o) get_spatial(o, id = "BOGUS"),
  get_colocalization  = function(o) get_colocalization(o, id = "BOGUS"),
  get_spatial_niche   = function(o) get_spatial_niche(o, id = "BOGUS"),
  get_rare_cells      = function(o) get_rare_cells(o, id = "BOGUS"),
  get_perturbation    = function(o) get_perturbation(o, id = "BOGUS"),
  get_detest          = function(o) get_detest(o, id = "BOGUS"),
  get_enrichment      = function(o) get_enrichment(o, id = "BOGUS"),
  get_mapping         = function(o) get_mapping(o, id = "BOGUS"),
  get_integration     = function(o) get_integration(o, id = "BOGUS"),
  get_annotation      = function(o) get_annotation(o, id = "BOGUS"),
  get_cci             = function(o) get_cci(o, id = "BOGUS")
)
for (nm in names(accs)) {
  r <- tryCatch(accs[[nm]](sce), error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1L && grepl("^ERROR", r)) {
    cat(sprintf("  %-22s %s\n", nm, r))
  } else if (is.null(r)) {
    cat(sprintf("  %-22s NULL (correct: no such id)\n", nm))
  } else {
    cat(sprintf("  %-22s *** WRONG RECORD RETURNED *** id=%s method=%s\n",
                nm, r$id %||% "(none)", r$method %||% "(none)"))
  }
}

cat("\n=========== E2  sclet_resolve_group_column() with no active identity ===========\n")
sce2 <- scuttle::mockSCE(ncells = 20, ngenes = 50)
cat("  ActiveIdent(sce) is NULL:", is.null(ActiveIdent(sce2)), "\n")
r2 <- tryCatch(sclet_resolve_group_column(sce2, as.data.frame(SummarizedExperiment::colData(sce2)),
                                           group = NULL),
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat("  sclet_resolve_group_column():", r2, "\n")
for (fn in c("RunCellChat", "RunCellPhoneDB")) {
  rr <- tryCatch(do.call(fn, list(sce2)), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  %-16s -> %s\n", fn, if (is.character(rr)) rr else "no error"))
}

cat("\n=========== E3  readiness-gate call graph ===========\n")
gates <- c("check_annotation_readiness", "check_cnv_readiness", "check_communication_readiness",
           "check_de_readiness", "check_enrichment_readiness", "check_fate_drivers_readiness",
           "check_fate_readiness", "check_integration_readiness", "check_milo_readiness",
           "check_multimodal_readiness", "check_rare_cell_readiness", "check_spatial_readiness",
           "check_trajectory_readiness", "check_velocity_readiness")
src_files <- list.files(file.path(pkg, "R"), pattern = "[.]R$", full.names = TRUE)
for (g in gates) {
  calls <- character(0)
  for (f in src_files) {
    ln <- readLines(f, warn = FALSE)
    hit <- grep(paste0("\\b", g, "\\s*\\("), ln)
    hit <- hit[!grepl("^\\s*#", ln[hit])]            # drop roxygen @examples
    hit <- hit[!grepl(paste0("<- function"), ln[hit])]  # drop the definition
    if (length(hit)) calls <- c(calls, sprintf("%s:%d", basename(f), hit))
  }
  cat(sprintf("  %-36s %s\n", g,
              if (length(calls)) paste(calls, collapse = ", ") else "*** NO CALL SITE ANYWHERE ***"))
}

cat("\n=========== E4  gate package check vs analysis package requirement ===========\n")
recs2 <- list()
for (f in src_files) {
  e <- parse(f, keep.source = TRUE)
  for (x in e) {
    if (!is.call(x) || !is.symbol(x[[1]])) next
    if (!(as.character(x[[1]]) %in% c("<-", "=") && length(x) >= 3)) next
    if (!is.symbol(x[[2]])) next
    b <- paste(deparse(x[[3]]), collapse = "\n")
    pk <- unique(sub('requireNamespace\\("', "", regmatches(b, gregexpr('requireNamespace\\("[^"]+"', b))))
    pk <- sub('"$', "", pk)
    recs2[[as.character(x[[2]])]] <- pk
  }
}
pairs <- list(
  c("check_de_readiness", "FindAllMarkers"),
  c("check_de_readiness", "FindMarkers"),
  c("check_velocity_readiness", "RunVelocity"),
  c("check_trajectory_readiness", "RunSlingshot"),
  c("check_integration_readiness", "RunIntegration"),
  c("check_spatial_readiness", "RunSpatialDeconvolution"),
  c("check_spatial_readiness", "RunSpatialColocalization"),
  c("check_annotation_readiness", "RunSingleR"),
  c("check_annotation_readiness", "RunKNNPredict"),
  c("check_annotation_readiness", "RunSymphonyMapping"),
  c("check_rare_cell_readiness", "RunRareCellDetection"),
  c("check_fate_readiness", "RunCellRank"),
  c("check_milo_readiness", "RunMilo"),
  c("check_enrichment_readiness", "RunEnrichment"),
  c("check_cnv_readiness", "RunCNVInference"),
  c("check_communication_readiness", "RunCellChat")
)
for (p in pairs) {
  g <- recs2[[p[1]]]; a <- recs2[[p[2]]]
  cat(sprintf("  %-30s checks: %-34s | %-24s requires: %s\n", p[1],
              if (length(g)) paste(g, collapse = ",") else "(none)",
              p[2], if (length(a)) paste(a, collapse = ",") else "(none)"))
}

cat("\n=========== E5  RunCCI keep_object -> RunNicheNet ===========\n")
cat("  formals(RunNicheNet):", paste(names(formals(RunNicheNet)), collapse = ", "), "\n")
cat("  has '...':", "..." %in% names(formals(RunNicheNet)),
    " has 'keep_object':", "keep_object" %in% names(formals(RunNicheNet)), "\n")
cat("  formals(RunCCI):", paste(names(formals(RunCCI)), collapse = ", "), "\n")
cat("  RunCCI body forwards keep_object to:",
    grep("keep_object", deparse(body(RunCCI)), value = TRUE), "\n")
