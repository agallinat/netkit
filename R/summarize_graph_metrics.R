#' Summarize Topological Properties of a Graph
#'
#' Computes a comprehensive set of global topological metrics for an input graph,
#' including basic structure, connectivity, spectral properties, and complexity.
#' Supports both `igraph` objects and data frames representing edge lists.
#'
#' @param graph An `igraph` object or a data frame with columns `from` and `to` representing an edge list.
#'
#' @return A tibble with one row and multiple columns, each representing a graph-level metric.
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
#' @importFrom igraph is_igraph graph_from_data_frame as.undirected degree V E vertex_attr as_adjacency_matrix
#' @importFrom igraph components induced_subgraph edge_density diameter
#' @importFrom igraph mean_distance transitivity assortativity_degree
#' @importFrom igraph betweenness vcount ecount modularity cluster_louvain
#' @importFrom tibble tibble
#' @importFrom ineq Gini
#' @importFrom Matrix Diagonal
#'
#' @references
#' - Newman, M. E. J. (2010). *Networks: An Introduction*. Oxford University Press.
#' - Estrada, E. (2012). *The Structure of Complex Networks: Theory and Applications*. Oxford University Press.
#' - Latora, V., Nicosia, V., & Russo, G. (2017). *Complex Networks: Principles, Methods and Applications*. Cambridge University Press.
#' - Louvain modularity method: Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008). *Fast unfolding of communities in large networks*. J. Stat. Mech., 2008(10), P10008.
#'
#' @examples
#' \dontrun{
#' g <- igraph::sample_gnp(200, 0.05, directed = F)
#' summarize_graph_metrics(g)
#' }
#'
#' @export
summarize_graph_metrics <- function(graph) {

  # --- Validate input ---
  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
  } else if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be either an igraph object or a data.frame representing an edge list.")
  }

  directed <- is_directed(graph)

  if (directed) {
    graph <- as.undirected(graph, mode = "collapse")
  }

  comps <- components(graph)
  lcc <- induced_subgraph(graph, which(comps$membership == which.max(comps$csize)))
  deg <- degree(graph)
  single_nodes <- sum(deg == 0)
  nodes <- vcount(graph)

  # Laplacian spectrum (sparse)
  A <- as_adjacency_matrix(graph, sparse = TRUE)
  L <- Diagonal(x = rowSums(A)) - A

  eigen_vals <- tryCatch({
    vals <- RSpectra::eigs(L, k = 2, which = "SM")$values
    sort(Re(vals))[2]
  }, error = function(e) NA)

  Algebraic_connectivity = eigen_vals

  # Entropy (degree)
  p <- table(deg) / length(deg)
  degree_entropy <- -sum(p * log2(p), na.rm = TRUE)

  # Gini coefficient of degree
  gini_deg <- ineq::Gini(deg)

  # Modularity (via Louvain, fast for large graphs)
  if (nodes > 2) {
    mod_score <- tryCatch({
      igraph::modularity(cluster_louvain(graph))
    }, error = function(e) NA)
  } else {
    mod_score <- NA
  }

  # Approximate betweenness
  Avg_betweenness <- if (nodes > 5000) NA else mean(betweenness(graph))

  # Approximate mean distance
  Average_path_length <- tryCatch({
    mean_distance(lcc)
  }, error = function(e) NA)

  data.frame(
    Nodes = nodes,
    Edges = ecount(graph),
    Is_directed = directed,
    Density = edge_density(graph),
    Diameter = diameter(lcc),
    Average_path_length = Average_path_length,
    Clustering_coefficient = transitivity(graph, type = "average"),
    Degree_assortativity = assortativity_degree(graph),
    Avg_degree = mean(deg),
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
