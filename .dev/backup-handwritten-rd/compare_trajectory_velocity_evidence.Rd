\name{compare_trajectory_velocity_evidence}
\alias{compare_trajectory_velocity_evidence}
\title{Compare trajectory and velocity pseudotime evidence summaries}
\usage{
compare_trajectory_velocity_evidence(object, trajectory_id = NULL, velocity_id = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}
\item{trajectory_id}{Optional analysis-run name for the trajectory evidence node (matching the \code{name} parameter of \code{run_trajectory}). When \code{NULL}, the single available trajectory run is used; if more than one exists, the function returns \code{status = "ambiguous_run"}.}
\item{velocity_id}{Optional analysis-run name for the velocity evidence node (matching the \code{name} parameter of \code{run_velocity}). When \code{NULL}, the single available velocity run is used; if more than one exists, the function returns \code{status = "ambiguous_run"}.}
}
\value{
A bounded list with \code{status}, spread comparisons, and explicit caveats.
When either evidence node is missing or incomplete, returns
\code{status = "not_available"} with a typed \code{reason}. When multiple runs
of either type exist and no id was given, returns \code{status = "ambiguous_run"}
listing the available ids in \code{available_trajectory_runs} and
\code{available_velocity_runs}. The result contains no per-cell values and
always sets \code{raw_values_included = FALSE}.
}
\description{
Reads only the bounded evidence nodes emitted by completed trajectory and
velocity runs (\code{ev:trajectory_<name>} and \code{ev:velocity_<name>}).
Compares their pseudotime spread summaries (relative IQR / max ratios) and
reports whether the two estimates are of similar magnitude. Never reads raw
per-cell \code{colData} columns, never claims either pseudotime is correct,
and never writes a new evidence node. The \code{claim_level} is always
\code{"consistent_with"}, matching both source evidence nodes.
}
