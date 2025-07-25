% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_attributes.R
\name{assign_attributes}
\alias{assign_attributes}
\title{Assign Vertex and Edge Attributes to an igraph Graph}
\usage{
assign_attributes(
  graph,
  nodes_table = NULL,
  edge_table = NULL,
  overwrite = TRUE
)
}
\arguments{
\item{graph}{An \code{igraph} object or a data frame containing a symbolic edge list.}

\item{nodes_table}{Optional. A \code{data.frame} whose first column corresponds to vertex names.}

\item{edge_table}{Optional. A \code{data.frame} whose first two columns correspond to source and target vertices.}

\item{overwrite}{Logical. If \code{TRUE}, existing attributes are overwritten. If \code{FALSE}, existing attributes are preserved. Default is \code{TRUE}.}
}
\value{
An \code{igraph} object with added or updated attributes.
}
\description{
Adds or updates vertex and edge attributes in an \code{igraph} object
using user-provided metadata tables. Vertex attributes are matched by the first column of \code{nodes_table},
and edge attributes are matched using the first two columns of \code{edge_table}, taking graph direction into account.
Only matching nodes and edges are updated. Warnings are issued when there are unmatched entries.
}
