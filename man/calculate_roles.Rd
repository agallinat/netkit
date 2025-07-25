% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_roles.R
\name{calculate_roles}
\alias{calculate_roles}
\title{Calculate Network Roles Based on Within-Module Z-Score and Participation Coefficient}
\usage{
calculate_roles(
  graph,
  communities = NULL,
  cluster.method = "spinglass",
  plot = TRUE,
  highlight_roles = TRUE,
  hub_z = 2.5,
  label_region = NULL,
  label.size = 12
)
}
\arguments{
\item{graph}{An \code{igraph} object representing the network.}

\item{communities}{Optional. A community clustering object (as returned by an \code{igraph} clustering function), or a named membership vector. If \code{NULL}, community detection is performed using \code{cluster.method}.}

\item{cluster.method}{Character. Clustering algorithm to use if \code{communities} is \code{NULL}. Default is \code{"spinglass"}. Passed to \code{find_modules()}.}

\item{plot}{Logical. Whether to generate a 2D plot of participation coefficient (P) vs. within-module z-score (z). Default is \code{TRUE}.}

\item{highlight_roles}{Logical. If \code{TRUE}, the role regions in the z–P plane are shaded for visual clarity. Default is \code{TRUE}.}

\item{hub_z}{Numeric. Threshold for defining hubs in terms of within-module z-score. Default is \code{2.5}.}

\item{label_region}{Optional character vector of role labels (e.g., \code{c("R4", "R7")}) indicating which role regions should have their nodes labeled in the plot. Default is \code{NULL}.}

\item{label.size}{Numeric. Base font size for plot text. Default is \code{12}.}
}
\value{
A list with three elements:
\describe{
\item{\code{plot}}{The \code{ggplot2} object (only if \code{plot = TRUE}).}
\item{\code{roles_definitions}}{A data frame describing the seven role types and their conditions.}
\item{\code{result}}{A data frame with node-level information: node name, module, z-score, participation coefficient, and assigned role.}
}
}
\description{
Implements the node role classification system of Guimerà & Amaral (2005) by calculating
the within-module degree z-score and the participation coefficient for each node in a network.
Nodes are assigned to one of seven role categories (R1–R7) based on their local modular connectivity.
}
\details{
If no community structure is provided, modules are automatically detected using the specified clustering method.
The function can optionally produce a 2D role plot (z vs. P) highlighting the canonical role regions.

The node roles are defined as:

\tabular{ll}{
R1 \tab Ultra-peripheral (non-hub): \eqn{z < 2.5, P <= 0.05} \cr
R2 \tab Peripheral (non-hub): \eqn{z < 2.5, 0.05 < P <= 0.6} \cr
R3 \tab Non-hub connector: \eqn{z < 2.5, 0.6 < P <= 0.8} \cr
R4 \tab Non-hub kinless: \eqn{z < 2.5, P > 0.8} \cr
R5 \tab Provincial hub: \eqn{z >= 2.5, P <= 0.3} \cr
R6 \tab Connector hub: \eqn{z >= 2.5, 0.3 < P <= 0.75} \cr
R7 \tab Kinless hub: \eqn{z >= 2.5, P > 0.75} \cr
}
}
\examples{
\dontrun{
g <- igraph::sample_gnp(200, 0.05, directed = F)
igraph::V(g)$name <- as.character(1:200)
result <- calculate_roles(g, plot = TRUE)
head(result$result)
}

}
\references{
Guimerà, R., & Amaral, L. A. N. (2005). Functional cartography of complex metabolic networks. \emph{Nature}, 433(7028), 895–900. \doi{10.1038/nature03288}
}
\seealso{
\code{\link[=find_modules]{find_modules()}}, \code{\link[igraph:cluster_spinglass]{igraph::cluster_spinglass()}}, \code{\link[igraph:communities]{igraph::membership()}}
}
