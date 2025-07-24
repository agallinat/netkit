#' Perform Network Diffusion from Seed Nodes
#'
#' Applies network diffusion techniques to propagate influence from a set of seed nodes across a graph.
#' Supports Laplacian smoothing, heat diffusion, and random walk with restart (RWR).
#'
#' @param graph An \code{igraph} object representing the network to analyze or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes. Must have named vertices.
#' @param seed_nodes Character vector of seed node names (must match `V(graph)$name`).
#' @param method Character. Diffusion method to use:
#'   * `"laplacian"`: Solves the linear system \eqn{(I + \alpha L)^{-1} f_0},
#'     where \eqn{L} is the (normalized) graph Laplacian and \eqn{\alpha} is a smoothing parameter.
#'     Internally, a sparse Cholesky decomposition is used for efficiency.
#'   * `"heat"`: Applies the heat diffusion model \eqn{e^{-tL} f_0}, where \eqn{t} controls diffusion time.
#'     If available, a sparse approximation method is used to avoid dense matrix exponential.
#'   * `"rwr"`: Random Walk with Restart. Iteratively solves \eqn{f = (1 - r)Pf + r f_0},
#'     where \eqn{P} is the transition matrix and \eqn{r} is the restart probability.
#' @param alpha Damping factor for the Laplacian method. Default is `0.7`.
#' @param t Time parameter for the heat diffusion method. Default is `1`.
#' @param restart_prob Restart probability (usually between 0.3 and 0.7) for the RWR method. Default is `0.3`.
#' @param normalize Logical. Whether to normalize the adjacency matrix (symmetric normalization for undirected graphs). Default is `TRUE`.
#' @param n_permutations Integer. Number of permutations to run for empirical p-value estimation (default `1000`).
#' @param seed Optional integer for reproducible random number generation. If `NULL` (default), seed is not set.
#' @param verbose Logical. If `TRUE` (default), displays a progress bar during permutations.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{node}{Node name}
#'   \item{score}{Diffusion score representing influence from the seed nodes}
#' }
#'
#' @details
#' This function allows flexible application of network diffusion strategies, useful in systems biology
#' (e.g., gene prioritization, pathway propagation), network analysis, and disease gene discovery. The underlying
#' matrix operations are based on well-established diffusion models from graph theory.
#'
#' For `"rwr"` (random walk with restart), the algorithm iteratively propagates scores until convergence
#' based on a row-normalized transition matrix. Recommended method for large networks.
#'
#' For `"laplacian"` and `"heat"`, the graph Laplacian is computed from the (optionally normalized) adjacency matrix.
#'
#' @references
#' Köhler S, Bauer S, Horn D, Robinson PN. Walking the interactome for prioritization of candidate disease genes.
#' \emph{Am J Hum Genet}. 2008;82(4):949–958. \doi{10.1016/j.ajhg.2008.02.013}
#'
#' Vanunu O, Magger O, Ruppin E, Shlomi T, Sharan R. Associating genes and protein complexes with disease via network propagation.
#' \emph{PLoS Comput Biol}. 2010;6(1):e1000641. \doi{10.1371/journal.pcbi.1000641}
#'
#' @examples
#' \dontrun{
#' g <- sample_gnp(100, 0.05, directed = F)
#' V(g)$name <- as.character(seq_len(vcount(g)))
#' seed_nodes <- sample(V(g)$name, 5)
#' network_diffusion_with_pvalues(g, seed_npodes, method = "laplacian")
#' }
#'
#' @importFrom igraph is_igraph V is_directed as_adjacency_matrix vertex_attr vertex_attr<- vcount
#' @importFrom Matrix Diagonal
#' @importFrom expm expm
#' @importFrom dplyr arrange desc
#'
#' @export
network_diffusion_with_pvalues <- function(graph,
                                           seed_nodes,
                                           method = c("laplacian", "heat", "rwr"),
                                           alpha = 0.7, t = 1, restart_prob = 0.3,
                                           normalize = TRUE,
                                           n_permutations = 1000,
                                           seed = NULL,
                                           verbose = TRUE) {

  set.seed(seed)
  all_nodes <- igraph::vertex_attr(graph, "name")
  n_seeds <- length(seed_nodes)
  seed_nodes <- intersect(seed_nodes, all_nodes)

  if (length(seed_nodes) == 0) stop("None of the seed nodes are in the graph.")

  # 1. Compute real diffusion scores
  real_scores_df <- network_diffusion(graph, seed_nodes, method, alpha, t, restart_prob, normalize)
  real_scores <- setNames(real_scores_df$score, real_scores_df$node)

  if (verbose) {
    message(sprintf("Running %d permutations with %d random seed nodes each (parallelized)...",
                    n_permutations, n_seeds))
  }

  # Precompute diffusion matrix
  precomp <- prepare_diffusion(graph = graph,
                               method = method,
                               alpha = alpha, t = t, restart_prob = restart_prob,
                               normalize = normalize)

  # 2. Prepare parallel plan
  future::plan("multisession")

  # 3. Run permutations in parallel
  perm_results <- future.apply::future_lapply(seq_len(n_permutations), function(i) {
    perm_seeds <- sample(setdiff(all_nodes, seed_nodes), n_seeds)
    null_df <- network_diffusion(graph, perm_seeds, method, alpha, t, restart_prob, normalize, precompute = precomp)
    setNames(null_df$score, null_df$node)
  }, future.seed = seed)

  # 4. Convert list of named vectors to a matrix
  null_scores_mat <- do.call(cbind, perm_results)
  null_scores_mat <- null_scores_mat[names(real_scores), , drop = FALSE]  # align rows

  # 5. Compute empirical p-values
  p_values <- mapply(function(real, null_dist) {
    (sum(null_dist >= real) + 1) / (length(null_dist) + 1)
  }, real = real_scores, null_dist = split(null_scores_mat, row(null_scores_mat)))

  # 6. Return result
  result <- tibble::tibble(
    node = names(real_scores),
    score = as.numeric(real_scores),
    p_empirical = as.numeric(p_values),
    stringsAsFactors = FALSE
  )
  result <- result[order(result$p_empirical), ]
  return(result)
}
