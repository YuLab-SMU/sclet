\name{summarize_multimodal_evidence}
\alias{summarize_multimodal_evidence}
\title{Summarize the registered multimodal modalities}
\usage{
summarize_multimodal_evidence(object, name = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{name}{Optional multimodal record id. When supplied, only that record is
summarized; when \code{NULL}, every multimodal record is summarized.}
}
\value{
A read-only summary list. When no modality is available and no record exists,
\code{status = "not_available"} with \code{reason = "no_multimodal_record"}.

Otherwise \code{status = "available"} with \code{n_modalities},
\code{modality_kinds} (plain assay-platform labels such as \code{"ADT"} and
\code{"ATAC"}), \code{modalities} (per-kind \code{n_features} and \code{n_cells}
as stored in the state record artifacts, each with
\code{raw_values_included = FALSE}), \code{has_ADT} and \code{has_ATAC} as the
workflow record stores them, \code{workflow_record_available},
\code{claim_level = "consistent_with"}, \code{caveats},
\code{recommendation = NULL} and \code{raw_values_included = FALSE}.

No cell barcode, antibody name, feature or gene symbol, file path, sample
identifier or per-cell vector is ever returned.
}
\description{
Returns a bounded, read-only summary of the registered multimodal state
records. Only counts, booleans, plain assay-platform labels and aggregate
statements are emitted.

The \code{caveats} state that modality registration records what was
\strong{attached} and tracked and performs no integration or batch correction;
that the counts describe the stored altExp dimensions rather than biological
content; and that no per-cell value or antibody/gene identity is included.

The function does not mutate the object, writes no evidence node, and never
recommends an integration method.
}
