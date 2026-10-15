\name{compare_multimodal_evidence}
\alias{compare_multimodal_evidence}
\title{Compare registered multimodal modalities}
\usage{
compare_multimodal_evidence(object)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}
}
\value{
A read-only comparison list. Fewer than two registered modalities return
\code{status = "not_available"} with
\code{reason = "at_least_two_modalities_required"} and \code{n_modalities}
reporting how many were found.

Otherwise \code{status = "available"} with \code{n_modalities},
\code{modality_codes} (anonymized codes matching \code{group_N}, derived from a
fixed ADT-then-ATAC ordering so the same modality always yields the same code),
\code{modalities} (per-code \code{n_features} and \code{n_cells}, each with
\code{raw_values_included = FALSE}), \code{claim_level = "consistent_with"},
\code{caveats}, \code{recommendation = NULL} and
\code{raw_values_included = FALSE}.
}
\description{
Observationally compares the registered multimodal modalities. The function
reports how many modalities are registered and how many features and cells each
one stores, and nothing else.

The \code{caveats} state that this compares the stored \strong{dimensions} of
different assay platforms rather than biological signal, that a larger modality
is not a better one, and that no modality is preferred, ranked or recommended.

The function never selects, ranks or prefers a modality, never describes one as
the informative or dominant modality, writes no evidence node, and removes or
filters no cell. Because the modalities are different assay platforms with
different feature spaces, the counts are not comparable as measures of
biological content.
}
