\name{summarize_milo_evidence}
\alias{summarize_milo_evidence}
\title{Summarize Milo results in bounded form without claiming a conclusion}
\usage{
summarize_milo_evidence(object, id = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{id}{Optional Milo record id. When \code{NULL}, the active Milo record is
used.}
}
\value{
A read-only summary list. When no Milo record exists, \code{status} is
\code{"not_available"} with a \code{reason} of \code{"no_milo_record"} or
\code{"no_milo_results"}. Otherwise \code{status} is \code{"available"} with
\code{n_da_results}, \code{n_significant}, FDR quantiles
(\code{fdr_quantile_25}, \code{fdr_median}, \code{fdr_quantile_75}),
\code{p_value_median}, \code{n_samples}, \code{n_conditions},
\code{n_contrasts}, \code{available_milo_ids},
\code{evidence_node_recorded}, \code{claim_level = "associated"},
\code{recommendation = NULL}, \code{raw_values_included = FALSE} and
\code{caveats}. A count that the record does not carry is reported as \code{NA},
not as zero.
}
\description{
Reports the same neighbourhood counts and FDR quantiles the Milo evidence node
carries, and writes nothing.
}
\details{
The function never ranks neighbourhoods, never names a sample or condition
label, and never recommends a design column, a contrast or a graph parameter.
Differential abundance results are associations between neighbourhoods and
conditions under a chosen graph, normalisation and test; they are neither causal
nor measured quantities, and a neighbourhood that did not reach significance is
not evidence of no difference, because the graph and the FDR procedure both
depend on the design the user declared. Only counts and quantiles are returned;
the DA table, the design frame and the model formula stay in the state record.
}
\seealso{\code{\link{check_milo_readiness}}, \code{\link{RunMilo}},
\code{\link{ConfirmAIDesignSemantics}}}
