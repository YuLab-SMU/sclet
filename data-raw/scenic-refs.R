## Download SCENIC reference files for bookdown CI.
## These are small text files (~1-2 MB total) that can run in CI.
## The cisTarget feather DB (~297 MB for human 10kb) is NOT downloaded;
## see bookdown/program-regulon-mechanism.Rmd for download instructions.
##
## Output:
##   ../data/allTFs_hg38.txt        -- TF list (text, ~15 KB)
##   ../data/motifs-v9-nr.hgnc.tbl  -- motif annotations (~1.5 MB)

tf_url <- "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"
motif_url <- "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

dir.create("../data", showWarnings = FALSE)

message("Downloading TF list...")
download.file(tf_url, "../data/allTFs_hg38.txt")

message("Downloading motif annotations...")
download.file(motif_url, "../data/motifs-v9-nr.hgnc.tbl")

message("Done. cisTarget feather DB (~297 MB) must be downloaded separately:")
message("  https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
