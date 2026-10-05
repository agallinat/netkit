#' Summarize Topological Properties of a Graph
#'
#' Computes a comprehensive set of global topological metrics for an input graph,
#' including basic structure, connectivity, spectral properties, and complexity.
#' Supports both `igraph` objects and data frames representing edge lists.
#'
#' @param graph An `igraph` object or a data frame with columns `from` and `to` representing an edge list.
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#'
#' @return A one-row `data.frame`, each column a graph-level metric.
#'
#'   Metrics that are mathematically undefined for the input are `NaN` rather than
#'   substituted values, as returned by the underlying \pkg{igraph} and \pkg{ineq}
#'   functions. On degenerate graphs this is expected: an edgeless graph has no paths
#'   (`Average_path_length`), no connected triples (`Clustering_coefficient`), no
#'   degree variance (`Degree_assortativity`) and a zero mean degree
#'   (`Gini_degree`).
#'
#' @details
#' Metrics computed:
#' \itemize{
#'   \item Number of nodes and edges
#'   \item Directed TRUE/FALSE
#'   \item Graph density
#'   \item Diameter and average path length of the largest connected component
#'   \item Clustering coefficient (transitivity)
#'   \item Degree assortativity
#'   \item Average degree and betweenness centrality
#'   \item Number of connected components and size of the largest connected component
#'   \item Number of single nodes
#'   \item Algebraic connectivity (second-smallest Laplacian eigenvalue)
#'   \item Degree entropy (Shannon entropy of the degree distribution)
#'   \item Gini coefficient of node degrees
#'   \item Modularity of the community structure (via Louvain algorithm)
#' }
#'
#' When `weights` are supplied, the metrics switch to their weighted definitions:
#' `Avg_degree` becomes mean vertex strength; `Diameter`, `Average_path_length`
#' and `Avg_betweenness` are computed on edge *costs*;
#' `Clustering_coefficient` uses Barrat's weighted transitivity;
#' `Degree_assortativity` correlates strengths rather than degrees;
#' `Degree_entropy` and `Gini_degree` describe the strength distribution; and the
#' Laplacian behind `Algebraic_connectivity` and the Louvain run behind
#' `Modularity` both use edge strengths. Two columns are added, `Is_weighted` and
#' `Avg_strength`, so a weighted and an unweighted summary can be row-bound.
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @importFrom igraph is_igraph graph_from_data_frame as_undirected degree V E vertex_attr as_adjacency_matrix
#' @importFrom igraph components induced_subgraph edge_density diameter
#' @importFrom igraph mean_distance transitivity assortativity_degree
#' @importFrom igraph betweenness vcount ecount modularity cluster_louvain
#' @importFrom igraph strength assortativity
#' @importFrom tibble tibble
#' @importFrom ineq Gini
#' @importFrom Matrix Diagonal
#'
#' @references
#' - Newman, M. E. J. (2010). *Networks: An Introduction*. Oxford University Press.
#' - Estrada, E. (2012). *The Structure of Complex Networks: Theory and Applications*. Oxford University Press.
#' - Latora, V., Nicosia, V., & Russo, G. (2017). *Complex Networks: Principles, Methods and Applications*. Cambridge University Press.
#' - Louvain modularity method: Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008). *Fast unfolding of communities in large networks*. J. Stat. Mech., 2008(10), P10008.
#' - Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A. (2004). *The architecture of complex weighted networks*. PNAS, 101(11), 3747-3752.
#'
#' @examples
#' g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
#' summarize_graph_metrics(g)
#'
#' # Weighted: declare what the attribute means, because igraph alone would read
#' # it as a cost for path metrics and as a strength for modularity.
#' igraph::E(g)$confidence <- runif(igraph::ecount(g), 0.1, 1)
#' summarize_graph_metrics(g, weights = "confidence", weight_type = "strength")
#'
#' @export
summarize_graph_metrics <- function(graph,
                                    weights = NULL,
                                    weight_type = c("strength", "distance")) {

  # --- Validate input ---
  graph <- as_netkit_graph(graph)
  w <- as_netkit_weights(graph, weights, weight_type)

  directed <- is_directed(graph)

  if (directed) {
    # Collapsing a directed graph merges reciprocal edges, so an edge-aligned
    # weight vector no longer lines up. Carry the weights on the graph itself and
    # let igraph sum them, then read them back off.
    if (w$weighted) {
      graph <- igraph::set_edge_attr(graph, ".netkit_s", value = w$strength)
      graph <- igraph::set_edge_attr(graph, ".netkit_d", value = w$distance)
      graph <- as_undirected(graph, mode = "collapse",
                             edge.attr.comb = list(.netkit_s = "sum",
                                                   .netkit_d = "min",
                                                   "ignore"))
      w$strength <- igraph::edge_attr(graph, ".netkit_s")
      w$distance <- igraph::edge_attr(graph, ".netkit_d")
    } else {
      graph <- as_undirected(graph, mode = "collapse")
    }
  }

  comps <- components(graph)
  lcc_ids <- which(comps$membership == which.max(comps$csize))
  lcc <- induced_subgraph(graph, lcc_ids)
  nodes <- vcount(graph)

  # Degree is a count; strength is the weighted analogue. Both are reported when
  # weights are in use, because the count is still the more interpretable of the
  # two and dropping it would make weighted and unweighted output incomparable.
  deg <- degree(graph)
  single_nodes <- sum(deg == 0)
  strengths <- if (w$weighted) {
    igraph::strength(graph, weights = w$strength)
  } else {
    deg
  }

  # Path-based measures take the *cost* vector, and they run on the largest
  # connected component rather than the whole graph. Rather than trying to work
  # out which positions of w$distance survive induced_subgraph(), carry the costs
  # as an edge attribute and let igraph subset them: that cannot silently
  # misalign, which a positional index very easily can.
  lcc_dist <- NULL
  if (w$weighted) {
    lcc <- induced_subgraph(
      igraph::set_edge_attr(graph, ".netkit_d", value = w$distance),
      lcc_ids
    )
    lcc_dist <- igraph::edge_attr(lcc, ".netkit_d")
  }

  # Laplacian spectrum (sparse). Built from strengths so that algebraic
  # connectivity measures how strongly, not merely how many ways, the graph holds
  # together.
  A <- if (w$weighted) {
    g_w <- igraph::set_edge_attr(graph, ".netkit_s", value = w$strength)
    igraph::as_adjacency_matrix(g_w, sparse = TRUE, attr = ".netkit_s")
  } else {
    as_adjacency_matrix(graph, sparse = TRUE)
  }
  L <- Diagonal(x = rowSums(A)) - A

  eigen_vals <- tryCatch({
    vals <- RSpectra::eigs(L, k = 2, which = "SM")$values
    sort(Re(vals))[2]
  }, error = function(e) NA)

  Algebraic_connectivity = eigen_vals

  # Entropy and Gini of the degree (or strength) distribution.
  p <- table(strengths) / length(strengths)
  degree_entropy <- -sum(p * log2(p), na.rm = TRUE)
  gini_deg <- ineq::Gini(strengths)

  # Modularity (via Louvain, fast for large graphs). Community detection reads a
  # weight as a strength, which is the one place igraph's implicit default
  # already agrees with netkit's.
  if (nodes > 2) {
    mod_score <- tryCatch({
      igraph::modularity(cluster_louvain(graph, weights = w$strength))
    }, error = function(e) NA)
  } else {
    mod_score <- NA
  }

  # Approximate betweenness
  Avg_betweenness <- if (nodes > 5000) {
    NA
  } else {
    mean(betweenness(graph, weights = w$distance))
  }

  # Approximate mean distance
  Average_path_length <- tryCatch({
    mean_distance(lcc, weights = lcc_dist)
  }, error = function(e) NA)

  Diameter <- tryCatch({
    diameter(lcc, weights = lcc_dist)
  }, error = function(e) NA)

  # Barrat's definition is the weighted generalization of transitivity; igraph
  # exposes it as type = "barrat" and it needs a strength per edge.
  Clustering_coefficient <- if (w$weighted) {
    tryCatch({
      mean(transitivity(graph, type = "barrat", weights = w$strength),
           na.rm = TRUE)
    }, error = function(e) NA)
  } else {
    transitivity(graph, type = "average")
  }

  # assortativity_degree() has no weights argument; the weighted analogue is the
  # same correlation computed over strengths.
  Degree_assortativity <- if (w$weighted) {
    tryCatch(igraph::assortativity(graph, values = strengths),
             error = function(e) NA)
  } else {
    assortativity_degree(graph)
  }

  data.frame(
    Nodes = nodes,
    Edges = ecount(graph),
    Is_directed = directed,
    Is_weighted = w$weighted,
    Density = edge_density(graph),
    Diameter = Diameter,
    Average_path_length = Average_path_length,
    Clustering_coefficient = Clustering_coefficient,
    Degree_assortativity = Degree_assortativity,
    Avg_degree = mean(deg),
    Avg_strength = mean(strengths),
    Avg_betweenness = Avg_betweenness,
    Components = comps$no,
    Single_nodes = single_nodes,
    LCC_size = max(comps$csize),
    LCC_percent = max(comps$csize) / nodes,
    Algebraic_connectivity = Algebraic_connectivity,
    Degree_entropy = degree_entropy,
    Gini_degree = gini_deg,
    Modularity = mod_score
  )
}
