## Build all cached bookdown data in one local entry point.
## Output files are written to ../data.
##
## This script is intended as the payload for a manual or scheduled cache
## refresh workflow; bookdown itself only consumes the saved objects.

scripts <- c(
    "pbmc-3k-sce.R",
    "pbmc-scvi.R",
    "hermann-sperm-velo.R",
    "hermann-sperm-dynamics.R",
    "visium-dlpfc-spatial.R",
    "scenic-refs.R"
)

for (script in scripts) {
    message("=== Running ", script, " ===")
    source(script, local = new.env(parent = globalenv()))
}
