\name{summarize_communication_evidence}
\alias{summarize_communication_evidence}
\title{Summarize CCI results in bounded form without claiming a conclusion}
\usage{
summarize_communication_evidence(object, id = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{id}{Optional communication record id. When \code{NULL}, the active
communication record is used.}
}
\value{
A read-only summary list. When no communication record or interaction table
exists, \code{status} is \code{"not_available"} with a \code{reason} of
\code{"no_communication_record"} or \code{"no_interaction_table"}. Otherwise
\code{status} is \code{"available"} with group-flow aggregates
(\code{n_interactions}, \code{n_source_groups}, \code{n_target_groups},
anonymized \code{source_group_codes} and \code{target_group_codes} aligned by
position with a row-major \code{pair_flow_counts_rowmajor} vector), strength
quantiles, \code{n_significant_pairs}, an engine \code{method_code},
\code{evidence_node_recorded}, \code{claim_level = "consistent_with"},
\code{recommendation = NULL}, \code{raw_values_included = FALSE} and
\code{caveats}.
}
\description{
Reports the same bounded group-flow aggregates the communication evidence node
carries, and writes nothing.
}
\details{
The function never ranks populations or pathways, never names a ligand or
receptor, and never recommends an engine. A predicted ligand-receptor pair is
not evidence that the cells signalled, and the absence of a pair is uninformative
because each engine filters against its own curated database and complex model.
Only group codes and numeric aggregates are returned; the full interaction
table, which carries gene identity, stays in the analysis layer and is read in R
through \code{\link{get_cci}}.
}
\seealso{\code{\link{check_communication_readiness}}, \code{\link{RunCCI}}}
