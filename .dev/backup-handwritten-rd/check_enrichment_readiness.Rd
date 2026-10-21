\name{check_enrichment_readiness}
\alias{check_enrichment_readiness}
\title{Check whether an object is ready for enrichment analysis}
\usage{
check_enrichment_readiness(object, detest_id = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{detest_id}{Optional detest record id. When supplied, that exact record
must be present.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks} (detest
record availability, resolved detest id, whether the marker table carries a
cluster column, cluster identity resolvability, resolved cluster key, package
availability, the accepted \code{available_dbs} values, and
\code{raw_values_included = FALSE}), \code{blocked_actions}
(\code{"enrichment"} when not ready), \code{questions} (human-readable
descriptions of missing prerequisites) and \code{notes} (an explicit statement
that enrichment categories are derived from the differential-expression marker
table and are not independent measurements).
}
\description{
Reports whether a completed differential-expression test and a usable cluster
assignment are present for \code{RunEnrichment()}. No database, ontology, key
type or significance cutoff is proposed, ranked or recommended: each is a
modelling choice reserved for the user.
}
\details{
The diagnostic also guards a silent failure specific to this stage. When the
marker table carries no cluster column, \code{RunEnrichment()} labels every
marker as a single pooled group, emits a warning rather than an error, and
returns a comparison of one set against itself. That outcome is reported here
as a blocking gate, because a meaningless result is worse for a planning agent
than a blocked one. No gene symbol and no enrichment category is read into the
return value.
}
