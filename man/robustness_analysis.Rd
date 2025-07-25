% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustness_analysis.R
\name{robustness_analysis}
\alias{robustness_analysis}
\title{Network Robustness Analysis via Node Removal Simulation}
\usage{
robustness_analysis(
  graph,
  removal_strategy = c("random", "degree", "betweenness"),
  steps = 50,
  metrics = c("lcc_size", "efficiency", "n_components"),
  n_reps = 50,
  plot = TRUE,
  seed = NULL
)
}
\arguments{
\item{graph}{An \code{igraph} object representing the network to plot or a
data frame containing a symbolic edge list in the first two columns. Additional
columns are considered as edge attributes. Must be undirected; directed graphs
will be converted.}

\item{removal_strategy}{Character. Strategy used for node removal. Options are:
\code{"random"}, \code{"degree"}, \code{"betweenness"}, or the name of a numeric vertex attribute.
Custom attributes are interpreted as priority scores (higher = removed first).}

\item{steps}{Integer. Number of removal steps (default: 50).}

\item{metrics}{Character vector. Structural metrics to compute at each step.
Options include: \code{"lcc_size"}, \code{"efficiency"}, and \code{"n_components"}.}

\item{n_reps}{Integer. Number of simulation repetitions (only relevant if \code{removal_strategy = "random"}).}

\item{plot}{Logical. If \code{TRUE}, a robustness plot is generated.}

\item{seed}{Integer or NULL. Random seed for reproducibility.}
}
\value{
A list with:
\describe{
\item{\code{all_results}}{A data frame with simulation results across all steps and repetitions.}
\item{\code{summary}}{A summarized data frame (mean and SD) if \code{n_reps > 1}, otherwise raw results.}
\item{\code{auc}}{Named list of AUC (area under the curve) values for each selected metric.}
}

If \code{plot = TRUE}, a ggplot object is generated showing the evolution of selected metrics
as nodes are progressively removed.
}
\description{
Simulates the removal of nodes from a network using various strategies and
evaluates how the structure degrades using selected robustness metrics.
Useful for assessing the vulnerability or resilience of a graph.
}
\details{
This function builds on classic approaches in network science for evaluating structural robustness
(e.g., Albert, Jeong, & Barabási, Nature 2000) by simulating progressive node removal and quantifying
the degradation of key topological features.

For deterministic strategies (\code{"degree"}, \code{"betweenness"}, or custom attributes),
nodes are removed in a fixed priority order. For the \code{"random"} strategy, the process is
repeated \code{n_reps} times, and the results are aggregated.

#' @references
Albert R, Jeong H, Barabási AL. Error and attack tolerance of complex networks.
Nature. 2000;406(6794):378–382. \doi{10.1038/35019019}

The function uses:
\itemize{
\item \strong{Largest Connected Component (LCC)}: Size of the largest remaining component.
\item \strong{Global Efficiency}: Average inverse shortest path length among all pairs.
\item \strong{Number of Components}: Total number of disconnected components.
}

Additionally, Area Under the Curve (AUC) is calculated for each metric, providing a scalar summary
of robustness. A higher AUC indicates greater resilience (i.e., slower degradation).

The implementation is inspired by principles described in:
Albert R, Jeong H, Barabási AL. \emph{Error and attack tolerance of complex networks}. Nature. 2000;406:378–382.
(\doi{10.1038/35019019})

The function uses:
\itemize{
\item Size of the largest connected component (\code{lcc_size})
\item Global efficiency (average inverse shortest path length)
\item Number of components (\code{n_components})
}
to evaluate how robust the network remains during progressive node failure.
}
\examples{
\dontrun{
g <- igraph::sample_pa(100)
robustness_analysis(g, removal_strategy = "degree", metrics = c("lcc_size", "n_components"))
}

}
