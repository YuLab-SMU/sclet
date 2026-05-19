#' @importFrom basilisk BasiliskEnvironment
sclet_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_env",
    pkgname = "sclet",
    packages = c(
        "scanpy==1.10.1",
        "anndata==0.10.7",
        "numpy==1.26.4",
        "pandas==2.2.2",
        "pyscenic==0.12.1",
        "loompy==3.0.7",
        "scvelo==0.3.2",
        "scvi-tools",
        "cellrank",
        "cell2location",
        "celloracle"
    )
)
