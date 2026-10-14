\name{RunIntegrationRoutes}
\alias{RunIntegrationRoutes}
\title{Run multiple integration routes and register comparable route metrics}
\usage{
RunIntegrationRoutes(
  object,
  design,
  routes = c("raw", "fastMNN"),
  confirm = c("ask", "always", "never"),
  baseline = "raw",
  ...
)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{design}{A named list describing the experimental design. At minimum it
must contain \code{batch}, naming the batch column in \code{colData(object)}.
It may also contain \code{condition} for biological-preservation summaries.}

\item{routes}{Character vector of integration routes to compare. Supported
values are \code{"raw"}, \code{"fastMNN"}, \code{"Harmony"}, and
\code{"scVI"}. The \code{"raw"} route is a recorded baseline and does not run
an integration method.}

\item{confirm}{Confirmation policy. \code{"ask"} returns a dry-run preview in
non-interactive sessions, \code{"always"} executes immediately, and
\code{"never"} cancels without modifying the object.}

\item{baseline}{Route id used as the comparison baseline. Defaults to
\code{"raw"}.}

\item{...}{Reserved for future extensions. Additional arguments are currently
ignored.}
}
\value{
A list with components including:
\itemize{
  \item \code{object}: the updated \code{SingleCellExperiment} when execution
  occurred, otherwise the input object unchanged.
  \item \code{plan}: a structured route-execution summary.
  \item \code{validation}: design and baseline validation fields.
  \item \code{preview}: a dry-run preview with route ids and estimated costs.
  \item \code{execution}: a list with \code{allowed}, \code{performed}, and
  \code{n_routes}.
  \item \code{status}: one of \code{"clarification_required"},
  \code{"cancelled"}, \code{"dry_run"}, or \code{"completed"}.
  \item \code{report}: per-route metric summaries and trade-off annotations.
}
}
\description{
\code{RunIntegrationRoutes()} is the execution facade for controlled,
AI-guided integration comparison. It can register an uncorrected baseline,
execute one or more supported integration methods, and write bounded,
comparable summaries back to the analysis ledger.
}
\details{
Each executed route is recorded as an integration analysis/state entry with a
uniform metric contract. The current metrics layer records route summaries for
batch mixing, biological preservation, cluster stability, and runtime. When a
metric cannot be computed, the route keeps an explicit
\code{not_available}/reason record instead of fabricating a value.

The function never auto-selects a best route. Any detected trade-off is reported
in \code{report$tradeoffs_detected}, while \code{report$recommendation} remains
\code{NULL}. This keeps route comparison observational rather than prescriptive.

If \code{design} is missing, not a list, or names a non-existent batch column,
the function returns \code{status = "clarification_required"} and does not run
integration.
}
\seealso{
\code{\link{RunIntegration}} for a single integration run, and
\code{\link{CompareAIAnalyses}} for read-only comparison of already recorded
route summaries.
}
