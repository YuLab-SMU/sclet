\name{advanced-ai-diagnostics}
\alias{summarize_qc_by_group}
\alias{summarize_pca_metadata_association}
\alias{summarize_cluster_sample_composition}
\alias{summarize_small_clusters}
\alias{check_integration_readiness}
\alias{check_annotation_readiness}
\alias{check_rare_cell_readiness}
\alias{summarize_small_cluster_evidence}
\alias{check_trajectory_readiness}
\alias{summarize_trajectory_cluster_order}
\title{Deterministic diagnostics for advanced AI analysis}
\usage{
summarize_qc_by_group(object, group)

summarize_pca_metadata_association(object, metadata, reduction = NULL)

summarize_cluster_sample_composition(object, sample, cluster = NULL)

summarize_small_clusters(object, threshold = 10L, cluster = NULL)

check_integration_readiness(object, design = NULL)

check_annotation_readiness(object, design = NULL)

check_rare_cell_readiness(object, cluster = NULL)

summarize_small_cluster_evidence(object, cluster = NULL, size_threshold = 10L,
  label = NULL)

check_trajectory_readiness(object, reduction = NULL)

summarize_trajectory_cluster_order(object, cluster = NULL, reduction = NULL,
  n_dims = 2L)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}
\item{group}{A \code{colData} column name or one value per cell.}
\item{metadata}{A \code{colData} column name or one value per cell.}
\item{reduction}{Optional reduction name; PCA is selected when available for
PCA association, and UMAP when available for trajectory checks.}
\item{sample}{A sample \code{colData} column name or one value per cell.}
\item{cluster}{Optional cluster \code{colData} column name or one value per cell.}
\item{threshold}{Maximum number of cells for a small cluster.}
\item{design}{Optional named list with sample, batch, condition, subject,
reference, labels, and group column names used by readiness checks.}
\item{size_threshold}{Maximum number of cells for a population to be
considered small.}
\item{label}{Optional single cluster label to summarize.}
\item{n_dims}{Number of leading embedding dimensions to summarize.}
}
\value{
Bounded, read-only structured summaries. Missing prerequisites return a result
with \code{status = "not_available"} or \code{"not_ready"} and a reason. Raw
metadata values and assay values are not included.
}
\description{
These functions compute deterministic summaries used by the advanced AI
investigation layer. They summarize QC by group, PCA association with metadata,
cluster/sample composition, small clusters, integration readiness,
annotation/marker readiness, rare-cell readiness, and trajectory readiness
without executing a modifying analysis.

\code{check_trajectory_readiness()} never proposes a trajectory root: choosing
which cluster represents the origin is a biological assumption reserved for the
user, so the result reports only whether a cluster assignment and a usable
embedding are present.
\code{summarize_trajectory_cluster_order()} reports the distribution of clusters
across leading embedding dimensions without ranking clusters or interpreting
cluster order as time.

\code{summarize_small_cluster_evidence()} aggregates the already computed
signals bearing on a small population (QC columns, doublet calls, recorded
marker/annotation analyses, sample or batch composition) and reports how many
independent signal classes are available. It deliberately returns no verdict
about whether a population is real or noise.
}
