\name{RunBasicWorkflow}
\alias{RunBasicWorkflow}
\title{Run the standard basic single-cell workflow}
\usage{
RunBasicWorkflow(
  object,
  n_features = 2000,
  n_pcs = 30,
  cluster_resolution = 0.5,
  steps = c("normalize", "variable_features", "scale", "pca", "neighbors",
    "clusters", "umap"),
  verbose = TRUE
)
}
\arguments{
\item{object}{a SingleCellExperiment object.}

\item{n_features}{number of highly variable features requested by
\code{FindVariableFeatures()}.}

\item{n_pcs}{number of principal components requested by \code{RunPCA()} and the
number of leading dimensions used by \code{FindNeighbors()} and \code{RunUMAP()}.}

\item{cluster_resolution}{clustering resolution passed to \code{FindClusters()}.}

\item{steps}{Character vector naming the steps to run, in execution order.
Supported names are \code{"normalize"}, \code{"variable_features"}, \code{"scale"},
\code{"pca"}, \code{"neighbors"}, \code{"clusters"} and \code{"umap"}.}

\item{verbose}{Logical. Report progress with \code{cli::cli_alert_info()}, matching
the logging style used elsewhere in the package.}
}
\value{
the updated SingleCellExperiment object.
}
\description{
A beginner-oriented convenience wrapper that chains the standard sclet
preprocessing, dimensionality reduction, graph and clustering steps in the
conventional order. It is purely deterministic: no AI call is made and no
action registry is involved. This is the "do not think about the details"
entry point; expert users should call the individual functions directly.
}
\details{
Every step delegates to the corresponding existing sclet function, which is
responsible for its own state registration. This function only sequences the
calls and never re-implements any analysis logic.

\code{RunBasicWorkflow()} is a human-facing helper and is deliberately not
registered in \code{AIDefaultExecutionRegistry()}.
}
\examples{
sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = matrix(rpois(50 * 40, lambda = 5), nrow = 50L, ncol = 40L))
)
sce <- RunBasicWorkflow(sce, n_pcs = 5)
SingleCellExperiment::reducedDimNames(sce)
}