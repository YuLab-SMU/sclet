\name{check_claim_ceiling}
\alias{check_claim_ceiling}
\title{Check what the recorded evidence can carry, before a claim is made}
\usage{
check_claim_ceiling(object, claims = NULL, evidence_ids = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{claims}{Optional list of proposed claims, each a list with
\code{claim_level} (character) and \code{evidence_refs} (character vector).
\code{NULL} reports only the ceiling of the object.}

\item{evidence_ids}{Optional restriction of which evidence nodes to consider.
Defaults to every completed evidence node on the object.}
}
\value{
A read-only summary list with elements \code{status}
(\code{"read_only"}), \code{ceiling_status} (the underlying ceiling status),
\code{ceiling_claim_level} (the strongest level the evidence can carry, or
\code{NULL} when nothing is supported), \code{ceiling_reasons}, the number of
independent \code{n_independent_groups}, \code{n_supporting_nodes}, and
\code{raw_values_included = FALSE}. When \code{claims} is supplied, a
\code{verdicts} list is added with one entry per claim, each carrying
\code{claim_level} (derived), \code{requested_claim_level} (the original,
preserved), \code{ceiling_status}, \code{ceiling_reasons} and a \code{status}
of \code{"allowed"}, \code{"downgraded"}, \code{"overclaimed"} or
\code{"unsupported"}.
}
\description{
Reports what the recorded evidence can actually carry, and whether a proposed
claim fits underneath that ceiling. \code{RecordAIResult(audit_claims = TRUE)}
enforces the ceiling at the write gate; this is the read-only half of the same
contract, for the moments when no claim exists yet.

The ceiling is computed by grouping the supporting nodes into lines of
independent support. A single line is capped at \code{"associated"}, so a lone
node can never be upgraded into a stronger statement; two or more independent
lines are required to reach \code{"consistent_with"}, and the ceiling is
bounded by the weakest supporting node either way. \code{user_decision} nodes
are excluded, because a human decision cannot support an AI-generated claim.

Three rules judge a proposed claim. \code{"causal"} is never an acceptable AI
claim level and is reported \code{"overclaimed"}, never silently downgraded. A
requested level above the ceiling is reported \code{"downgraded"} with the
corrected level attached, so the downgrade is visible rather than silent. A
claim citing no resolvable evidence is reported \code{"unsupported"} and
derives to \code{"hypothesis"}, the only level that explicitly self-limits and
asserts nothing. The model's original request is preserved in
\code{requested_claim_level} in every case.

The ceiling is the strongest the evidence can carry: this function only ever
lowers or preserves a requested level and never raises one. It records no
evidence, writes no finding and modifies nothing on the object.
}
