% Please edit documentation in R/ai-integration-routes.R
\name{RunIntegrationRoutes}
\alias{RunIntegrationRoutes}
\title{Run multiple integration strategies and register comparable evidence metrics}
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

\item{design}{A named list describing the experiment design with at least a
\code{batch} element naming the batch column in \code{colData}. May optionally
include a \code{condition} element naming the biological condition column.}

\item{routes}{Character vector of integration strategies to run. Supported
routes are \code{"raw"} (baseline, no correction), \code{"fastMNN"},
\code{"Harmony"}, and \code{"scVI"}. Defaults to \code{c("raw", "fastMNN")}.}

\item{confirm}{Confirmation gate. \code{"ask"} performs a dry run in
non-interactive sessions; \code{"always"} executes; \code{"never"} returns a
cancelled result without modifying the object.}

\item{baseline}{Route id used as the comparison baseline for metrics.
Defaults to \code{"raw"}.}

\item{...}{Additional design clarification questions returned alongside
\code{status = "clarification_required"} when the design is incomplete.}
}
\value{
A list with the following elements:
\itemize{
  \item \code{object}: The updated \code{SingleCellExperiment} with each route
    recorded as a state record of type \code{integration}.
  \item \code{plan}: A structured plan summary listing routes and baseline.
  \item \code{validation}: Design validation fields.
  \item \code{preview}: Dry-run preview summary including estimated costs.
  \item \code{execution}: A list with \code{allowed}, \code{performed}, and
    \code{n_routes}.
  \item \code{status}: One of \code{"clarification_required"},
    \code{"cancelled"}, \code{"dry_run"}, or \code{"completed"}.
  \item \code{report}: Per-route metric summaries plus
    \code{tradeoffs_detected} and a \code{recommendation} that is always
    \code{NULL} in the current implementation (routes are never
    auto-selected).
}
}
\description{
\code{RunIntegrationRoutes} is the execution entry point for the AI-guided
multi-route integration phase. It registers a baseline route, executes the
requested integration strategies via the existing \code{RunIntegration}
methods, and records a uniform metric summary (batch mixing, biological
preservation, cluster stability, runtime) for each route so that
\code{CompareAIAnalyses()} can compare them without side effects. When the
design is missing or the batch column is absent the function returns
\code{status = "clarification_required"}.
}
