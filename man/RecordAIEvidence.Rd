\name{RecordAIEvidence}
\alias{RecordAIEvidence}
\alias{ValidateAIEvidenceRefs}
\alias{sclet_ai_evidence_independence}
\title{Record and validate bounded AI evidence nodes}
\usage{
RecordAIEvidence(object, evidence, source = NULL, parents = NULL, scope = NULL)

ValidateAIEvidenceRefs(
  object,
  refs,
  fingerprint = NULL,
  requesting_scope = NULL,
  requesting_kind = NULL
)

sclet_ai_evidence_independence(object, evidence_ids)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{evidence}{A structured evidence node with a unique \code{id} and supported \code{kind}.}

\item{source}{Optional completed analysis id supporting the evidence.}

\item{parents}{Optional evidence ids or dependency ids used to derive the
deterministic \code{dependency_group} hash together with \code{scope}.}

\item{scope}{Optional aggregate scope, such as anonymized groups.}

\item{refs}{Character evidence ids to validate.}

\item{fingerprint}{Optional source fingerprint to require.}

\item{requesting_scope}{Optional aggregate scope of the referencing context.
Evidence references must not scope outside the supporting evidence's own scope.}

\item{requesting_kind}{Optional kind of the referencing context, used to
enforce cross-kind rules such as refusing \code{user_decision} citations from
AI-generated claims.}

\item{evidence_ids}{Character vector of exactly two evidence ids to compare.}
}
\value{
\code{RecordAIEvidence()} returns the SCE with an \code{ai_evidence} state
record (plus a deterministic \code{dependency_group} derived from
\code{parents} + \code{scope}). \code{ValidateAIEvidenceRefs()} returns a
structured result with \code{valid}, \code{resolved}, and \code{errors}
fields plus \code{requesting_scope} and \code{requesting_kind} audit
metadata. \code{sclet_ai_evidence_independence()} returns a list with
\code{independent}, \code{shared_dependency_groups}, and
\code{lineage_paths}.
}
\description{
These functions provide a bounded evidence contract for advanced AI analysis.
Evidence values must be aggregate (private identifiers and complete matrices
are rejected). Each evidence node carries a deterministic
\code{dependency_group} derived from its \code{parents} and \code{scope} so
downstream comparisons can detect shared lineage. References are validated
against the current analysis ledger and must additionally satisfy: (a) the
referencing scope is contained in each evidence node's own scope; and (b)
AI-generated claims (\code{deterministic_summary}, \code{plot}, \code{test},
\code{state}) are not allowed to cite \code{user_decision} nodes directly.
}
