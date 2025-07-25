% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_hubs.R
\name{find_hubs}
\alias{find_hubs}
\title{Identify Hub Nodes in a Network}
\usage{
find_hubs(
  graph,
  method = c("zscore", "quantile"),
  degree_threshold = 3,
  betweenness_threshold = 1,
  degree_quantile = 0.95,
  betweenness_quantile = 0.95,
  log_transform = TRUE,
  plot = TRUE,
  focus_color = "darkgreen",
  label.size = 12,
  hub_names = TRUE,
  hub_cex = 3,
  gg_extra = list()
)
}
\arguments{
\item{graph}{An \code{igraph} object or a data frame containing a symbolic edge list in the
first two columns. Additional columns are considered as edge attributes.}

\item{method}{Character. Method to identify hubs: \code{"zscore"} (standardized metrics)
or \code{"quantile"} (empirical percentiles).}

\item{degree_threshold}{Numeric. Threshold for standardized degree (only used if \code{method = "zscore"}).}

\item{betweenness_threshold}{Numeric. Threshold for standardized betweenness (only used if \code{method = "zscore"}).}

\item{degree_quantile}{Numeric between 0 and 1. Quantile threshold for degree (used if \code{method = "quantile"}).}

\item{betweenness_quantile}{Numeric between 0 and 1. Quantile threshold for betweenness (used if \code{method = "quantile"}).}

\item{log_transform}{Logical. If \code{TRUE}, applies log-transformation to degree and betweenness metrics.}

\item{plot}{Logical. If \code{TRUE}, generates a plot of the degree vs. betweenness classification.}

\item{focus_color}{Character. Color to display in the focus area of the plot (hubs region).}

\item{label.size}{Numeric. Base font size for plot elements. Passed to \code{theme_classic}.}

\item{hub_names}{Logical. If \code{TRUE}, adds node labels to identified hubs on the plot.}

\item{hub_cex}{Numeric. Font size scaling factor for hub labels on the plot.}

\item{gg_extra}{List. Additional user-defined layers for the returned ggplot.
eg. list(ylim(-2,2), theme_bw(), theme(legend.position = "none"))}
}
\value{
A list with the following components:
\describe{
\item{\code{method}}{Description of the method and thresholds used.}
\item{\code{result}}{A \code{tibble} with node name, degree, betweenness, transformed metrics, and hub status.}
\item{\code{graph}}{The original graph with a new vertex attribute \code{is_hub}.}
}
If \code{plot = TRUE}, a scatter plot of degree vs. betweenness is displayed with
hub nodes highlighted.
}
\description{
This function identifies hub nodes in an \code{igraph} network based on
degree and betweenness centrality using either z-score or quantile thresholds.
Optionally, it visualizes the classification using a scatter plot with marginal histograms.
}
\examples{
\dontrun{
library(igraph)
g <- sample_pa(100)
find_hubs(g, method = "quantile", plot = TRUE)
}

}
