% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_diffusion.R
\name{prepare_diffusion}
\alias{prepare_diffusion}
\title{Prepare Diffusion Matrix}
\usage{
prepare_diffusion(
  graph,
  method = c("laplacian", "heat", "rwr"),
  alpha = 0.7,
  t = 1,
  restart_prob = 0.3,
  normalize = TRUE
)
}
\arguments{
\item{graph}{An \code{igraph} object or a data frame containing a symbolic edge list in the
first two columns. Additional columns are considered as edge attributes.}

\item{method}{Character string: one of \code{"laplacian"}, \code{"heat"}, or \code{"rwr"}.}

\item{alpha}{Numeric (used in \code{"laplacian"}).}

\item{t}{Time parameter (used in \code{"heat"}).}

\item{restart_prob}{Restart probability (used in \code{"rwr"}).}
}
\value{
A matrix representing the diffusion kernel.
}
\description{
Prepares and normalizes the diffusion kernel matrix to be used in network diffusion.
}
\keyword{internal}
