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
#'   A kernel built with different edge weights than requested is rejected rather
#'   than reused.
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. Signal propagates along edge *strengths*. See
#'   [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#' @param seed_weights Optional numeric vector of initial values for the seeds,
#'   replacing the default binary indicator. Either named (matched to
#'   `seed_nodes` by name) or unnamed and parallel to `seed_nodes`. Use this to
#'   diffuse from a continuous signal -- log fold changes, scores, prior
#'   probabilities -- rather than from set membership, which is what the
#'   propagation literature generally assumes.
#'
#' @return A tibble with two columns, sorted by descending score:
#' \describe{
#'   \item{node}{Node name}
#'   \item{score}{Diffusion score representing influence from the seed nodes}
#' }
#'
#' @inheritSection netkit-weights Edge weights
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
#' g <- igraph::sample_gnp(80, 0.06, directed = FALSE)
#' igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
#' seed_nodes <- igraph::V(g)$name[1:5]
#'
#' network_diffusion(g, seed_nodes, method = "laplacian")
#'
#' # Reuse a precomputed kernel across repeated calls.
#' kernel <- prepare_diffusion(g, method = "rwr")
#' network_diffusion(g, seed_nodes, method = "rwr", precompute = kernel)
#'
#' # Diffuse along edge strengths rather than treating every edge alike.
#' igraph::E(g)$confidence <- runif(igraph::ecount(g), 0.1, 1)
#' network_diffusion(g, seed_nodes, method = "rwr", weights = "confidence")
#'
#' # Start from a continuous signal instead of set membership.
#' network_diffusion(g, seed_nodes, method = "rwr",
#'                   seed_weights = c(2.4, -1.1, 0.7, 3.0, 0.2))
#'
#' @importFrom igraph is_igraph V is_directed as_adjacency_matrix vertex_attr vertex_attr<- vcount
#' @importFrom Matrix Diagonal Cholesky rowSums solve
#' @importFrom dplyr arrange desc
#'
#' @export
network_diffusion <- function(graph, seed_nodes,
                              method = c("laplacian", "heat", "rwr"),
                              alpha = 0.7, t = 1, restart_prob = 0.3,
                              normalize = TRUE, precompute = NULL,
                              weights = NULL,
                              weight_type = c("strength", "distance"),
                              seed_weights = NULL) {

  method <- match.arg(method)
  weight_type <- match.arg(weight_type)

  # Convert graph to igraph if needed
  graph <- as_netkit_graph(graph, backfill_names = TRUE)

  all_nodes <- vertex_attr(graph, "name")
  n <- length(all_nodes)

  # Seed vector. Binary by default; `seed_weights` lets the caller start from a
  # continuous signal instead (log fold changes, scores, prior probabilities),
  # which is what the propagation literature assumes and what diffusion from a
  # differential-expression result actually needs.
  f0 <- make_seed_vector(all_nodes, seed_nodes, seed_weights)

  # If precompute is provided, use it
  if (!is.null(precompute)) {
    if (precompute$method != method) stop("Requested diffusion method do not match pre-computed data.")

    # A kernel built from different weights describes a different graph. Without
    # this check, passing a weighted kernel and then forgetting `weights =` would
    # silently diffuse over the wrong matrix and still return plausible scores.
    requested_key <- weights_key(as_netkit_weights(graph, weights, weight_type,
                                                   warn_unused = FALSE))
    if (!is.null(precompute$weights_key) && precompute$weights_key != requested_key) {
      stop("Pre-computed kernel was built with different edge weights (",
           precompute$weights_key, ") than requested (", requested_key,
           "). Rebuild it with prepare_diffusion(), or pass the same ",
           "'weights' and 'weight_type' used to build it.", call. = FALSE)
    }

    L <- precompute$L
    ch <- precompute$ch
    P <- precompute$P
    use_sparse_P <- precompute$use_sparse_P
  } else {
    # Delegate to prepare_diffusion() rather than rebuilding the kernel inline.
    # The two used to be separate copies of the adjacency-normalize -> Laplacian
    # -> Cholesky/transition-matrix block, which meant any change to the
    # diffusion maths had to be made twice.
    kernel <- prepare_diffusion(
      graph = graph,
      method = method,
      alpha = alpha,
      t = t,
      restart_prob = restart_prob,
      normalize = normalize,
      weights = weights,
      weight_type = weight_type
    )
    L <- kernel$L
    ch <- kernel$ch
    P <- kernel$P
    use_sparse_P <- kernel$use_sparse_P
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
                # The convergence test is relative to the magnitude of the seed
                # vector, not absolute. With a binary f0 the two are equivalent,
                # but `seed_weights` lets the caller start from values of any
                # scale (log fold changes, read counts), and an absolute
                # threshold would then mean the achieved precision depended on
                # the units of the input.
                tol <- 1e-6 * max(1, sum(abs(f0)))
                max_iter <- 1000L

                f_prev <- f0
                converged <- FALSE
                for (iter in seq_len(max_iter)) {
                  f_new <- (1 - restart_prob) * (P %*% f_prev) + restart_prob * f0
                  if (sum(abs(f_new - f_prev), na.rm = TRUE) < tol) {
                    converged <- TRUE
                    break
                  }
                  f_prev <- f_new
                }
                # The loop used to be an unbounded `repeat`, which cannot
                # terminate if the iteration fails to contract.
                if (!converged) {
                  warning(sprintf(
                    "Random walk with restart did not converge in %d iterations.",
                    max_iter
                  ), call. = FALSE)
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
#' @param normalize Logical. Whether to symmetrically normalize the adjacency
#'   matrix. Default is `TRUE`.
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. Diffusion propagates along edge *strengths*, so a
#'   `"distance"` type is converted before use. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#'
#' @return A list with the kernel components: `L` (the Laplacian), `ch` (its
#'   Cholesky factor, for `"laplacian"`), `P` (the transition matrix, for
#'   `"rwr"`), `use_sparse_P`, `method`, and `weights_key` -- a digest of the
#'   weights the kernel was built from, which [network_diffusion()] checks
#'   before reusing it.
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @examples
#' g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
#' igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
#'
#' kernel <- prepare_diffusion(g, method = "laplacian")
#' kernel$method
#'
#' # Passing the kernel back in skips rebuilding it on every call.
#' network_diffusion(g, seed_nodes = c("1", "2"), method = "laplacian",
#'                   precompute = kernel)
#'
#' # A weighted kernel propagates along edge strengths.
#' igraph::E(g)$confidence <- runif(igraph::ecount(g), 0.1, 1)
#' w_kernel <- prepare_diffusion(g, method = "rwr", weights = "confidence")
#'
#' @keywords internal
#'
#' @export
#'
prepare_diffusion <- function(graph,
                              method = c("laplacian", "heat", "rwr"),
                              alpha = 0.7, t = 1, restart_prob = 0.3,
                              normalize = TRUE,
                              weights = NULL,
                              weight_type = c("strength", "distance")) {

  method <- match.arg(method)

  graph <- as_netkit_graph(graph, backfill_names = TRUE)
  w <- as_netkit_weights(graph, weights, weight_type)

  # Diffusion propagates along edge strengths: a stronger edge carries more
  # signal. Building the adjacency matrix without `attr` discarded weights
  # entirely, which made every weighted diffusion bit-identical to the
  # unweighted one -- silently, since the scores still looked reasonable.
  A <- if (w$weighted) {
    g_w <- igraph::set_edge_attr(graph, ".netkit_strength", value = w$strength)
    igraph::as_adjacency_matrix(g_w, sparse = TRUE, attr = ".netkit_strength")
  } else {
    as_adjacency_matrix(graph, sparse = TRUE)
  }
  is_directed_graph <- is_directed(graph)

  if (normalize && !is_directed_graph) {
    deg <- Matrix::rowSums(A)
    deg[deg == 0] <- 1
    D_inv_sqrt <- Diagonal(x = 1 / sqrt(deg))
    A <- D_inv_sqrt %*% A %*% D_inv_sqrt
  } else if (normalize && is_directed_graph) {
    # Previously only network_diffusion() warned here. Now that it delegates the
    # kernel construction to this function, the warning has to live here too.
    warning("Normalization for directed graphs is not supported. Skipping normalization.")
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

  # The kernel is reusable only for the weights it was built from. Carrying a key
  # lets network_diffusion() reject a stale one rather than silently diffusing
  # over the wrong matrix -- the same guard the `method` check already provides.
  list(L = L, ch = ch, P = P, use_sparse_P = use_sparse_P, method = method,
       weights_key = weights_key(w))
}

#' Summarize a resolved weight spec for kernel-reuse checks
#'
#' Internal. Two kernels are interchangeable only if they were built from the
#' same weights. Hashing the strength vector (rather than storing it) keeps the
#' kernel small while still detecting a mismatch.
#'
#' @param w The list returned by `as_netkit_weights()`.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
weights_key <- function(w) {
  if (!w$weighted) {
    return("unweighted")
  }
  paste0("weighted:", w$type, ":", length(w$strength), ":",
         format(sum(w$strength), digits = 15), ":",
         format(sum(w$strength^2), digits = 15))
}
