\name{ResolveAIClarifications}
\alias{ResolveAIClarifications}
\title{Resolve AI clarification questions interactively}
\usage{
ResolveAIClarifications(
  object,
  clarification,
  ask = NULL,
  apply_answer = NULL,
  record = TRUE
)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{clarification}{A clarification list as returned by
\code{sclet_ai_format_clarification()}, or \code{result$report$clarification} from
\code{RunAIAnalysis()}.}

\item{ask}{Optional function taking a question list and returning the user's
answer as a character scalar. Defaults to a terminal prompt. Supply this to embed
a different UI or to script the interaction.}

\item{apply_answer}{Optional function \code{(object, question, answer)} returning
the updated object, used for question types that have no built-in confirmation
function.}

\item{record}{Logical. Record accepted answers as \code{user_decision} evidence.}
}
\value{
A list with the updated \code{object}, a \code{status}, and a \code{transcript} of
per-question \code{resolved} / \code{skipped} / \code{invalid} outcomes.
}
\description{
Presents the structured questions produced by \code{sclet_ai_format_clarification()}
to the user, collects an answer for each, applies the answer through the existing
confirmation mechanism where one exists, and records each accepted answer as a
\code{user_decision} evidence node. This closes the loop described in the
clarification contract: the AI was blocked, the user was asked, and the user
answered.
}
\details{
The function never answers on the user's behalf. It proposes no default column,
no default cluster and no default reference, and an empty answer is treated as
"not answered" rather than being filled in. Answers are only applied when the
question type has a real confirmation function (\code{design_batch} goes through
\code{ConfirmAIDesignSemantics()}); other question types describe plan parameters,
so they are recorded for the audit trail and returned for the caller to feed back
into its plan.

In a non-interactive session without an \code{ask} function the function refuses to
prompt and returns \code{status = "needs_interactive"} without modifying the object
or recording anything.
}
\examples{
sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = matrix(1, 4, 4))
)
SummarizedExperiment::colData(sce)$batch <- rep(c("a", "b"), 2)
clarification <- sclet:::sclet_ai_format_clarification(
    "prerequisites not met for action i: design_semantics_not_confirmed: x"
)
resolved <- ResolveAIClarifications(
    sce,
    clarification,
    ask = function(question) "batch"
)
resolved$status
}