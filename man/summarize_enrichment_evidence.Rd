\name{summarize_enrichment_evidence}
\alias{summarize_enrichment_evidence}
\title{Summarize recorded enrichment evidence without claiming a result}
\usage{
summarize_enrichment_evidence(object, name = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{name}{Optional enrichment record id. When \code{NULL}, the active
enrichment record is used.}
}
\value{
A read-only summary list. When no enrichment record exists, \code{status} is
\code{"not_available"} with a \code{reason} of \code{"no_enrichment_record"} or
\code{"no_enrichment_result"}. Otherwise \code{status} is \code{"available"}
with bounded counts (\code{n_enriched}, \code{n_clusters_enriched}),
anonymized \code{cluster_codes} aligned by position with
\code{enriched_per_cluster}, adjusted-p-value quantiles, the resolved
\code{detest_id} and its \code{lineage_parents}, \code{recommendation = NULL},
\code{claim_level = "consistent_with"}, \code{raw_values_included = FALSE} and
\code{caveats}.
}
\description{
Reports what an enrichment run recorded, in the same bounded shape its evidence
node uses, and writes nothing.
}
\details{
The function never ranks categories, never recommends a database or ontology,
and never surfaces gene or term identity: the enriched category names, their
accessions and the driving genes stay in the state record, reachable in R
through \code{\link{get_enrichment}}. Because enrichment re-describes a single
differential-expression measurement, its evidence declares that measurement as
its parent, so a claim ceiling counts the two as one independent line of
support rather than two.
}
