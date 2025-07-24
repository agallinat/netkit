#' Perform Network Diffusion from Seed Nodes
#'
#' Applies network diffusion techniques to propagate influence from a set of seed nodes across a graph.
#' Supports Laplacian smoothing, heat diffusion, and random walk with restart (RWR).
#'
#' @param graph An \code{igraph} object representing the network to analyze or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes. Must have named vertices.
#' @param seed_nodes Character vector of seed node names (must match \code{V(graph)$name}).
#' @param method Character. Diffusion method to use:
#'   \itemize{
#'     \item \code{"laplacian"}: Solves the linear system \eqn{(I + \alpha L)^{-1} f_0},
#'       where \eqn{L} is the (normalized) graph Laplacian and \eqn{\alpha} is a smoothing parameter.
#'       Internally, a sparse Cholesky decomposition is used for efficiency.
#'     \item \code{"heat"}: Applies the heat diffusion model \eqn{e^{-tL} f_0}, where \eqn{t} controls diffusion time.
#'       A truncated Taylor expansion is used for approximation.
#'     \item \code{"rwr"}: Random Walk with Restart. Iteratively solves \eqn{f = (1 - r)Pf + r f_0},
#'       where \eqn{P} is the transition matrix and \eqn{r} is the restart probability.
#'   }
#' @param alpha Damping factor for the Laplacian method. Default is \code{0.7}.
#' @param t Time parameter for the heat diffusion method. Default is \code{1}.
#' @param restart_prob Restart probability (usually between 0.3 and 0.7) for the RWR method. Default is \code{0.3}.
#' @param normalize Logical. Whether to normalize the adjacency matrix (symmetric normalization for undirected graphs). Default is \code{TRUE}.
#' @param precompute Optional list of precompute diffusion matrices (e.g., Laplacian, Cholesky factor, or transition matrix).
#'   Use \code{\link{prepare_diffusion}()} to generate this object and avoid redundant computations when calling this function repeatedly (e.g., in greedy optimization).
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{node}{Node name}
#'   \item{score}{Diffusion score representing influence from the seed nodes}
#' }
#'
#' @details
#' This function allows flexible application of network diffusion strategies, useful in systems biology
#' (e.g., gene prioritization, pathway propagation), network analysis, and disease gene discovery.
#' The underlying matrix operations are based on well-established diffusion models from graph theory.
#'
#' For \code{"rwr"} (random walk with restart), the algorithm iteratively propagates scores until convergence
#' based on a row-normalized transition matrix. Recommended method for large networks.
#'
#' For \code{"laplacian"} and \code{"heat"}, the graph Laplacian is computed from the (optionally normalized) adjacency matrix.
#' For efficiency in iterative applications, precompute Laplacian and Cholesky decomposition using \code{\link{prepare_diffusion}()}.
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
#' network_diffusion(g, seed_npodes, method = "laplacian")
#' }
#'
#' @importFrom igraph is_igraph V is_directed as_adjacency_matrix vertex_attr vertex_attr<- vcount
#' @importFrom Matrix Diagonal Cholesky rowSums solve
#' @importFrom dplyr arrange desc
#'
#' @export
network_diffusion <- function(graph, seed_nodes,
                              method = c("laplacian", "heat", "rwr"),
                              alpha = 0.7, t = 1, restart_prob = 0.3,
                              normalize = TRUE, precompute = NULL) {

  method <- match.arg(method)

  # Convert graph to igraph if needed
  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
  } else if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be either an igraph object or a data.frame representing an edge list.")
  }

  if (is.null(vertex_attr(graph, "name"))) {
    vertex_attr(graph, "name") <- as.character(seq_along(1:vcount(graph)))
  }

  all_nodes <- vertex_attr(graph, "name")
  n <- length(all_nodes)

  # Seed vector
  f0 <- as.numeric(all_nodes %in% seed_nodes)
  names(f0) <- all_nodes

  # If precompute is provided, use it
  if (!is.null(precompute)) {
    if (precompute$method != method) stop("Requested diffusion method do not match pre-computed data.")
    L <- precompute$L
    ch <- precompute$ch
    P <- precompute$P
    use_sparse_P <- precompute$use_sparse_P
  } else {
    # Adjacency matrix
    A <- as_adjacency_matrix(graph, sparse = TRUE)

    is_directed_graph <- is_directed(graph)

    if (normalize && !is_directed_graph) {
      deg <- Matrix::rowSums(A)
      deg[deg == 0] <- 1
      D_inv_sqrt <- Diagonal(x = 1 / sqrt(deg))
      A <- D_inv_sqrt %*% A %*% D_inv_sqrt
    } else if (normalize && is_directed_graph) {
      warning("Normalization for directed graphs is not supported. Skipping normalization.")
    }

    # Laplacian
    D <- Diagonal(x = Matrix::rowSums(A))
    L <- D - A

    # Cholesky for Laplacian diffusion
    ch <- NULL
    if (method == "laplacian") {
      I <- Diagonal(n = nrow(L))
      ch <- Matrix::Cholesky(I + alpha * L, LDL = FALSE, perm = TRUE)
    }

    # Transition matrix for RWR (sparse)
    P <- NULL
    use_sparse_P <- FALSE
    if (method == "rwr") {
      row_sums <- Matrix::rowSums(A)
      row_sums[row_sums == 0] <- 1
      D_inv <- Diagonal(x = 1 / row_sums)
      P <- D_inv %*% A  # sparse row-normalized transition matrix
      use_sparse_P <- TRUE
    }
  }

  # Heat diffusion: exp(-tL) * f0 via truncated Taylor
  approx_expmv <- function(L, f0, t = 1, K = 20) {
    result <- f0
    term <- f0
    for (k in 1:K) {
      term <- (-t / k) * (L %*% term)
      result <- result + term
    }
    result
  }

  # Diffuse
  f <- switch(method,
              "laplacian" = {
                if (is.null(ch)) stop("Missing Cholesky decomposition for Laplacian diffusion.")
                Matrix::solve(ch, f0)
              },
              "heat" = {
                approx_expmv(L, f0, t = t, K = 20)
              },
              "rwr" = {
                f_prev <- f0
                repeat {
                  f_new <- (1 - restart_prob) * (P %*% f_prev) + restart_prob * f0
                  if (sum(abs(f_new - f_prev), na.rm = TRUE) < 1e-6) break
                  f_prev <- f_new
                }
                f_new
              }
  )

  tibble::tibble(
    node = all_nodes,
    score = as.numeric(f)
  ) |> dplyr::arrange(desc(score))
}

#'
#' Prepare Diffusion Matrix
#'
#' Prepares and normalizes the diffusion kernel matrix to be used in network diffusion.
#'
#' @param graph An \code{igraph} object or a data frame containing a symbolic edge list in the
#'   first two columns. Additional columns are considered as edge attributes.
#' @param method Character string: one of `"laplacian"`, `"heat"`, or `"rwr"`.
#' @param alpha Numeric (used in `"laplacian"`).
#' @param t Time parameter (used in `"heat"`).
#' @param restart_prob Restart probability (used in `"rwr"`).
#'
#' @return A matrix representing the diffusion kernel.
#' @keywords internal
#'
#' @export
#'
prepare_diffusion <- function(graph,
                              method = c("laplacian", "heat", "rwr"),
                              alpha = 0.7, t = 1, restart_prob = 0.3,
                              normalize = TRUE) {

  method <- match.arg(method)

  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
  } else if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be either an igraph object or a data.frame representing an edge list.")
  }

  if (is.null(vertex_attr(graph, "name"))) {
    vertex_attr(graph, "name") <- as.character(seq_along(1:vcount(graph)))
  }

  A <- as_adjacency_matrix(graph, sparse = TRUE)
  is_directed_graph <- is_directed(graph)

  if (normalize && !is_directed_graph) {
    deg <- Matrix::rowSums(A)
    deg[deg == 0] <- 1
    D_inv_sqrt <- Diagonal(x = 1 / sqrt(deg))
    A <- D_inv_sqrt %*% A %*% D_inv_sqrt
  }

  D <- Diagonal(x = Matrix::rowSums(A))
  L <- D - A

  ch <- NULL
  if (method == "laplacian") {
    I <- Diagonal(n = nrow(L))
    ch <- Matrix::Cholesky(I + alpha * L, LDL = FALSE, perm = TRUE)
  }

  P <- NULL
  use_sparse_P <- FALSE
  if (method == "rwr") {
    row_sums <- Matrix::rowSums(A)
    row_sums[row_sums == 0] <- 1
    D_inv <- Diagonal(x = 1 / row_sums)
    P <- D_inv %*% A
    use_sparse_P <- TRUE
  }

  list(L = L, ch = ch, P = P, use_sparse_P = use_sparse_P, method = method)
}
