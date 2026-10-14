\name{compare_trajectory_roots}
\alias{compare_trajectory_roots}
\title{Compare completed trajectory runs across explicit root hypotheses}
\usage{
compare_trajectory_roots(object, ids = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}
\item{ids}{Optional character vector of completed trajectory analysis ids to compare. If \code{NULL}, all completed slingshot trajectory runs with bounded evidence are used.}
}
\value{
A bounded list with \code{status}, per-run summaries, root-hypothesis summaries,
pairwise comparisons, and explicit caveats. Fewer than two usable runs return
\code{status = "not_available"} with a typed \code{reason}. The result contains
no per-cell values, always sets \code{recommendation = NULL}, and always sets
\code{raw_values_included = FALSE}.
}
\description{
Reads only completed slingshot trajectory state records and their bounded
\code{ev:trajectory_<name>} evidence nodes. Compares multiple trajectory runs
generated from different explicit \code{start_cluster} assumptions, summarizes
how many runs support each anonymized root code, and reports pairwise
differences in lineage count and relative pseudotime spread. Never reads
per-cell pseudotime vectors, never chooses a preferred root, and never writes a
new evidence node.
}
