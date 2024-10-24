## code to load the databases for runSCENIC

# load the defaultDbNames.RData from SCENIC
load("./inst/ext/defaultDbNames.RData")
usethis::use_data(defaultDbNames, overwrite = TRUE)
 
# load the motifAnootation of three species from RcisTarget
load("./inst/ext/motifAnnotations_dmel_v9.RData")
motifAnnotations_dmel <- motifAnnotations_dmel_v9
load("./inst/ext/motifAnnotations_mgi_v9.RData")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
load("./inst/ext/motifAnnotations_hgnc_v9.RData")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
usethis::use_data(motifAnnotations_dmel, overwrite = TRUE)
usethis::use_data(motifAnnotations_mgi, overwrite = TRUE)
usethis::use_data(motifAnnotations_hgnc, overwrite = TRUE)

