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
#' @param null Which permutation null to compare against:
#'   \describe{
#'     \item{`"degree_matched"`}{The default. Each real seed is replaced by a
#'       randomly chosen vertex of similar degree, so the permuted sets have the
#'       same degree profile as the real one.}
#'     \item{`"uniform"`}{Permuted seed sets are drawn uniformly from the
#'       non-seed vertices, matching only the seed count.}
#'   }
#'   See Details for why the default is degree-matched.
#' @param match_pool Integer. For `null = "degree_matched"`, how many
#'   nearest-degree candidates each real seed may be replaced by. Default is
#'   `10`. Smaller values match degree more tightly but give the null less
#'   variance; `Inf` makes every candidate eligible and so reduces exactly to
#'   `null = "uniform"`.
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#' @param seed_weights Optional numeric vector of initial seed values, passed to
#'   [network_diffusion()], replacing the default binary indicator. Note that the
#'   permutation null re-uses these magnitudes on the permuted seed sets, so the
#'   null tests the *position* of the seeds rather than their values.
#'
#' @inheritSection netkit-weights Edge weights
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
#' @section Parallel execution:
#' The permutation null is evaluated with [future.apply::future_lapply()], which
#' runs under whichever \pkg{future} plan is currently active. This function does
#' not set a plan itself, so by default permutations run **sequentially**. To
#' parallelize, set a plan once in your own session before calling:
#'
#' ```
#' future::plan("multisession", workers = 4)
#' res <- network_diffusion_with_pvalues(g, seed_nodes, n_permutations = 1000)
#' future::plan("sequential")   # release the workers when finished
#' ```
#'
#' Permutations are reproducible regardless of the plan: `seed` is passed to
#' `future_lapply(future.seed = )`, which generates parallel-safe RNG streams.
#'
#' @references
#' Köhler S, Bauer S, Horn D, Robinson PN. Walking the interactome for prioritization of candidate disease genes.
#' \emph{Am J Hum Genet}. 2008;82(4):949–958. \doi{10.1016/j.ajhg.2008.02.013}
#'
#' Vanunu O, Magger O, Ruppin E, Shlomi T, Sharan R. Associating genes and protein complexes with disease via network propagation.
#' \emph{PLoS Comput Biol}. 2010;6(1):e1000641. \doi{10.1371/journal.pcbi.1000641}
#'
#' @examples
#' g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
#' igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
#' seed_nodes <- igraph::V(g)$name[1:5]
#'
#' # n_permutations is reduced from its default of 1000 to keep the example fast;
#' # use the default or higher for real analyses.
#' network_diffusion_with_pvalues(g, seed_nodes, method = "laplacian",
#'                                n_permutations = 50, seed = 1, verbose = FALSE)
#'
#' @section Choice of null:
#'
#' Diffusion scores are strongly degree-dependent: a high-degree vertex receives
#' signal from more directions, so it scores highly almost regardless of where
#' the seeds are. A null that draws seed sets uniformly therefore compares a real
#' seed set -- which in most applications is enriched for well-studied, highly
#' connected vertices -- against a null composed mostly of poorly connected ones,
#' and reports much of the degree difference as significance.
#'
#' `null = "degree_matched"` removes that confound by replacing each real seed
#' with a randomly chosen vertex drawn from the `match_pool` candidates closest
#' to it in degree. The resulting p-values answer "are these seeds unusually well
#' *placed*?" rather than "are these seeds unusually well *connected*?", which is
#' almost always the intended question.
#'
#' How closely the match can be made is limited by the graph. If the seed set is
#' extreme -- the six highest-degree vertices of a scale-free network, say --
#' there are simply no other vertices of comparable degree to permute with, and
#' no permutation null can fully correct for degree. The function warns when the
#' closest available candidates are still far from the seeds, so the residual
#' confound is visible rather than assumed away.
#'
#' @importFrom igraph is_igraph V is_directed as_adjacency_matrix vertex_attr vertex_attr<- vcount degree
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
                                           verbose = TRUE,
                                           weights = NULL,
                                           weight_type = c("strength", "distance"),
                                           seed_weights = NULL,
                                           null = c("degree_matched", "uniform"),
                                           match_pool = 10) {

  method <- match.arg(method)
  weight_type <- match.arg(weight_type)
  null <- match.arg(null)

  # Only seed when asked: set.seed(NULL) re-seeds from the clock, so the
  # documented default silently destroyed the caller's RNG stream. Same defect
  # as robustness_analysis() had.
  if (!is.null(seed)) set.seed(seed)

  # This function is documented to accept a data.frame edge list but was the one
  # export that never routed through the shared input validator -- it called
  # igraph::vertex_attr() on the raw argument, so an edge list failed with
  # igraph's own "Must provide a graph object" rather than working.
  graph <- as_netkit_graph(graph, backfill_names = TRUE)

  all_nodes <- igraph::vertex_attr(graph, "name")
  requested_seeds <- as.character(seed_nodes)
  seed_nodes <- intersect(requested_seeds, all_nodes)

  if (length(seed_nodes) == 0) stop("None of the seed nodes are in the graph.")

  # Count the seeds that were actually used, not the ones that were requested.
  # n_seeds used to be set before this intersection, so any seed absent from the
  # graph made every permuted set *larger* than the real one -- which inflates
  # the null scores and biases every p-value, without anything to show for it.
  n_seeds <- length(seed_nodes)

  if (length(seed_nodes) < length(requested_seeds)) {
    warning(sprintf(
      "%d of %d seed nodes are not vertices of the graph and were dropped.",
      length(requested_seeds) - length(seed_nodes), length(requested_seeds)
    ), call. = FALSE)
  }

  # Resolve seed_weights to a plain positional vector aligned with the retained
  # seeds, so the same magnitudes can be reused on each permuted set. A named
  # vector keyed by the real seed names could not be matched to permuted names.
  seed_weights <- resolve_perm_seed_weights(seed_weights, requested_seeds,
                                            seed_nodes)

  # 1. Compute real diffusion scores
  real_scores_df <- network_diffusion(graph, seed_nodes, method = method,
                                      alpha = alpha, t = t,
                                      restart_prob = restart_prob,
                                      normalize = normalize,
                                      weights = weights,
                                      weight_type = weight_type,
                                      seed_weights = seed_weights)
  real_scores <- setNames(real_scores_df$score, real_scores_df$node)

  if (verbose) {
    # class(plan()) is e.g. c("FutureStrategy", "sequential", "uniprocess", ...);
    # element 2 is the strategy name ("sequential", "multisession", ...).
    plan_name <- class(future::plan())[2]
    message(sprintf("Running %d permutations with %d random seed nodes each (future plan: %s)...",
                    n_permutations, n_seeds, plan_name))
  }

  # Precompute diffusion matrix
  precomp <- prepare_diffusion(graph = graph,
                               method = method,
                               alpha = alpha, t = t, restart_prob = restart_prob,
                               normalize = normalize,
                               weights = weights,
                               weight_type = weight_type)

  # 2. Run permutations under whatever future plan the caller has set.
  #
  # This deliberately does NOT call future::plan(). A package must not change the
  # user's plan: doing so overrides their choice globally, and the "multisession"
  # workers it starts are never shut down, leaving socket connections open. That
  # surfaces as "checking examples ... ERROR / connections left open" under
  # R CMD check. The default plan is sequential; see @details for opting in to
  # parallel execution.

  # 3. Run permutations
  # Candidate pool and, for the degree-matched null, the bin structure. Both are
  # computed once outside the permutation loop.
  candidates <- setdiff(all_nodes, seed_nodes)
  if (length(candidates) < n_seeds) {
    stop(sprintf(
      paste0("Only %d non-seed vertices are available but %d are needed per ",
             "permutation. The seed set is too large a share of the graph for ",
             "a permutation null."),
      length(candidates), n_seeds
    ), call. = FALSE)
  }

  sampler <- if (null == "degree_matched") {
    make_degree_matched_sampler(graph, seed_nodes, candidates, match_pool)
  } else {
    function() sample(candidates, n_seeds)
  }

  perm_results <- future.apply::future_lapply(seq_len(n_permutations), function(i) {
    perm_seeds <- sampler()
    null_df <- network_diffusion(graph, perm_seeds, method = method,
                                 alpha = alpha, t = t,
                                 restart_prob = restart_prob,
                                 normalize = normalize, precompute = precomp,
                                 weights = weights, weight_type = weight_type,
                                 seed_weights = seed_weights)
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
    p_empirical = as.numeric(p_values)
  )
  result <- result[order(result$p_empirical), ]
  return(result)
}

#' Build a degree-matched seed sampler
#'
#' Internal. Returns a closure that draws a permuted seed set whose degrees match
#' the real seed set's, by replacing each real seed with a randomly chosen
#' candidate of similar degree.
#'
#' Nearest-degree matching rather than binning. Two binning schemes were tried
#' first and both matched only loosely on exactly the graphs this correction
#' matters for. Quantile bins of degree *value* collapse on a heavy-tailed
#' distribution -- most vertices of a scale-free graph have degree 1, so most
#' quantiles coincide and the surviving top bin spans everything from degree 3 to
#' the largest hub. Quantile bins of degree *rank* are equal-sized by
#' construction, but still wide: ten bins over 200 vertices put twenty vertices
#' per bin, and in the top bin their degrees ranged from 1 to 17, so the bin mean
#' was nowhere near the seed it was meant to match. Matching each seed
#' individually removes the binning granularity from the problem entirely.
#'
#' @param graph An `igraph` object.
#' @param seed_nodes The real seed names.
#' @param candidates Vertex names eligible as permuted seeds.
#' @param match_pool How many nearest-degree candidates each real seed may be
#'   replaced by. Small values match degree tightly but give the null little
#'   variance; `Inf` makes every candidate eligible and so reduces exactly to the
#'   uniform null.
#'
#' @return A function of no arguments returning a character vector of seed names.
#'
#' @keywords internal
#' @noRd
make_degree_matched_sampler <- function(graph, seed_nodes, candidates,
                                        match_pool) {

  if (!is.numeric(match_pool) || length(match_pool) != 1 || match_pool < 1) {
    stop("'match_pool' must be a single number of at least 1.", call. = FALSE)
  }

  deg <- igraph::degree(graph)
  names(deg) <- igraph::V(graph)$name

  # When every candidate is eligible for every seed, this *is* the uniform null.
  # Returning that sampler rather than an equivalent-in-distribution one makes
  # the equivalence exact, which is what the documentation claims, and avoids
  # warning below about a match that was never requested.
  if (match_pool >= length(candidates)) {
    n_seeds <- length(seed_nodes)
    return(function() sample(candidates, n_seeds))
  }

  cand_deg <- deg[candidates]
  k <- min(match_pool, length(candidates))

  # The eligible pool for each seed, computed once rather than per permutation.
  pools <- lapply(seed_nodes, function(s) {
    ord <- order(abs(cand_deg - deg[[s]]))
    candidates[utils::head(ord, k)]
  })

  # How well the match can possibly do is limited by what the graph contains: an
  # extreme seed set has no comparable vertices to be replaced by. Say so once,
  # rather than leaving the caller to assume the correction was exact.
  best <- vapply(seq_along(seed_nodes), function(i) {
    mean(deg[pools[[i]]])
  }, numeric(1))
  seed_deg <- deg[seed_nodes]
  if (any(seed_deg > 0) &&
      mean(abs(best - seed_deg)) > 0.5 * mean(seed_deg)) {
    warning(sprintf(
      paste0("The graph contains few vertices of comparable degree to the seed ",
             "set (mean seed degree %.1f, mean degree of the closest available ",
             "candidates %.1f), so the null is only approximately ",
             "degree-matched. This is a property of the graph, not of the ",
             "sampling: an extreme seed set has no comparable vertices to be ",
             "permuted with."),
      mean(seed_deg), mean(best)
    ), call. = FALSE)
  }

  function() {
    drawn <- character(0)
    for (pl in pools) {
      available <- setdiff(pl, drawn)
      if (length(available) == 0) {
        # The pool is exhausted because an earlier seed took its members. Widen
        # to every remaining candidate rather than failing, and take the closest
        # in degree so the match degrades gracefully.
        others <- setdiff(candidates, drawn)
        if (length(others) == 0) {
          stop("Not enough candidate vertices to build a degree-matched ",
               "permutation. Use null = \"uniform\", or reduce the seed set.",
               call. = FALSE)
        }
        available <- others[order(abs(deg[others] - mean(deg[pl])))][1]
      }
      drawn <- c(drawn, available[sample.int(length(available), 1)])
    }
    drawn
  }
}
