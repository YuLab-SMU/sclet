#' @importFrom basilisk BasiliskEnvironment
sclet_scenic_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_scenic_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "pyscenic=0.12.1",
        "loompy=3.0.7"
    )
)

sclet_scvi_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_scvi_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "scvi-tools=1.4.2"
    )
)

sclet_cellrank_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_cellrank_env",
    pkgname = "sclet",
    packages = c(
        "python=3.12",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "scvelo=0.3.2",
        "cellrank=2.0.7"
    )
)

sclet_cell2location_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_cell2location_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "scvi-tools=1.2.2",
        "cell2location=0.1.5"
    )
)

sclet_celloracle_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_celloracle_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "celloracle=0.20.0"
    )
)
