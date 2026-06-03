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
old_timeout <- getOption("timeout")
options(timeout = 300)

download_safe <- function(url, dest) {
    tryCatch({
        message("Downloading ", basename(dest), "...")
        download.file(url, dest)
        message("  OK: ", dest)
    }, error = function(e) {
        message("  SKIP: ", conditionMessage(e))
        message("  (SCENIC refs are optional; SCENIC chunks use eval=FALSE)")
    })
}

download_safe(tf_url, "../data/allTFs_hg38.txt")
download_safe(motif_url, "../data/motifs-v9-nr.hgnc.tbl")

options(timeout = old_timeout)
