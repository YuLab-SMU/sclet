\name{summarize_spatial_evidence}
\alias{summarize_spatial_evidence}
\title{Summarize the stored spatial analyses}
\usage{
summarize_spatial_evidence(object, analysis = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{analysis}{Optional analysis key to summarize. One of
\code{"spatial_deconv"}, \code{"colocalization"} or
\code{"spatial_niche"}. When \code{NULL}, every analysis present is
summarized.}
}
\value{
A read-only summary list with elements \code{status} (\code{"available"} or
\code{"not_available"}), \code{reason} (\code{"no_spatial_record"} when no
spatial state record exists, \code{"analysis_not_available"} when the
requested analysis has no record), \code{n_analyses}, \code{analyses}
(per-analysis bounded counts keyed by anonymized \code{route_N} codes, each
carrying \code{analysis_code}, \code{n_records}, \code{n_features},
\code{n_columns}, \code{n_cells}, the algorithm \code{methods},
\code{claim_level = "consistent_with"} and
\code{raw_values_included = FALSE}), \code{claim_level},
\code{caveats} and \code{raw_values_included = FALSE}.
}
\description{
Reads the spatial state records and returns a bounded summary of what they
hold. Feature column names, per-spot model estimates and input identifiers
exist in the underlying records but never leave this function: only record
counts, algorithm names, anonymized analysis codes and fixed caveats are
returned. No cell/spot barcode, sample identifier, file path, gene symbol,
raw matrix or per-cell/per-spot vector is emitted, no spatial method is
ranked or recommended, and no AI evidence is recorded.
}
