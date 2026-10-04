\name{AIInvestigate}
\alias{AIInvestigate}
\title{Investigate an advanced single-cell question without executing actions}
\usage{
AIInvestigate(object, question, design = NULL,
    scope = c("diagnosis", "advanced", "report"), model = NULL,
    structured_output = TRUE, fallback_on_structure_error = TRUE, ...)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}
\item{question}{A single non-empty scientific question.}
\item{design}{Optional named list describing sample, batch, condition, or subject semantics.}
\item{scope}{Investigation scope.}
\item{model}{Optional aisdk model. If omitted, the configured default model is used.}
\item{structured_output}{Logical; request structured provider output.}
\item{fallback_on_structure_error}{Logical; retain a normalized answer if structured output is unavailable.}
\item{...}{Additional adapter arguments.}
}
\value{
A read-only object of class \code{sclet_ai_result}. It contains findings,
evidence, warnings, recommendations, and proposed routes. It never executes
an action or mutates the input object.
}
\description{
\code{AIInvestigate()} combines a bounded dataset profile, the analysis ledger,
and deterministic readiness diagnostics for advanced analysis questions. If
sample, batch, condition, or subject semantics are not confirmed, the result
reports a clarification requirement instead of executing a dependent action.
}
