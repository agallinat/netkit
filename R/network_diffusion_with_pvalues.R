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
#' @param match_pool For `null = "degree_matched"`, the smallest number of
#'   candidates any seed may be drawn from. Candidates are grouped into adaptive
#'   degree bins holding at least this many members, and each seed is replaced
#'   by a vertex drawn uniformly from its own bin.
#'
#'   `NULL` (the default) sets it to `10 * length(seed_nodes)`, capped at the
#'   number of candidates. The seed count has to enter the choice: every
#'   permutation draws `length(seed_nodes)` members out of a bin, so a fixed
#'   small pool means consecutive permutations are largely the same vertices and
#'   the null collapses. A number sets it explicitly and warns if it is small
#'   enough to put the null in that regime. `Inf` puts every candidate in one
#'   bin and so reduces exactly to `null = "uniform"`.
#'
#'   Tighter matching is not automatically better: see the "Choice of null"
#'   section for the trade-off it makes against the validity of the p-values.
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
#' @return A tibble with one row per vertex, ordered by `p_empirical` and then by
#'   descending `z`:
#' \describe{
#'   \item{node}{Node name.}
#'   \item{score}{Diffusion score representing influence from the seed nodes.}
#'   \item{p_empirical}{One-sided empirical p-value; see "Interpreting the
#'     p-values".}
#'   \item{null_mean, null_sd}{Mean and standard deviation of the vertex's own
#'     null distribution.}
#'   \item{z}{`(score - null_mean) / null_sd`, or `NA` where the null has no
#'     variance. This is the column to rank on when many p-values tie at the
#'     floor, which with the default `n_permutations` they routinely do.}
#'   \item{n_valid}{Number of permutations the vertex was actually tested
#'     against, which is `n_permutations` minus the permutations that seeded the
#'     vertex itself. The p-value floor is `1 / (n_valid + 1)`, so this is what
#'     sets the best attainable significance for that vertex.}
#' }
#'
#' @section Interpreting the p-values:
#'
#' `p_empirical` is **one-sided**: it counts null scores at least as *high* as
#' the observed one, and so tests whether a vertex receives more signal than
#' chance placement of the seeds would give. There is no two-sided variant,
#' because a vertex receiving unusually little diffusion signal is not a
#' meaningful finding -- most vertices receive almost none.
#'
#' It uses the `(r + 1) / (n + 1)` convention and is therefore never zero. With
#' `n_valid` permutations the smallest attainable value is `1 / (n_valid + 1)`,
#' so the default of 1000 permutations cannot distinguish the strongest vertex
#' in the graph from the hundredth strongest -- they all report 0.000999. Rank
#' on `z` when that happens; it is why the column is returned.
#'
#' A vertex is tested only against the permutations in which it was **not**
#' itself a permuted seed. A permutation that happened to seed vertex *v*
#' deposits a large score at *v* for exactly the reason the real seeds top the
#' real run, and says nothing about whether *v* is unusually well placed with
#' respect to the real seed set; including it inflated the upper tail of *v*'s
#' null and made its p-value conservative. The effect fell hardest on the
#' vertices drawn most often, so it was not a uniform shift. `n_valid` reports
#' what each vertex was actually compared against.
#'
#' The real seed nodes are a special case worth knowing about. They are never
#' drawn as permuted seeds, so they are compared against a null in which they
#' are never seeded, while in the real run they are. They will therefore sit at
#' the top of the ranking essentially by construction, and their p-values are
#' not evidence of anything. Drop them before interpreting the table.
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
#' @section Performance:
#'
#' The permutations are not diffused one at a time. They are grouped into
#' batches and each batch is solved as a single multiple-right-hand-side
#' problem, which for `"laplacian"` means one sparse Cholesky solve against a
#' matrix of seed vectors rather than one solve per permutation. On a
#' 6000-vertex, 78000-edge graph with 200 permutations this took the whole call
#' from 6.8 s to 3.4 s for `"laplacian"`, 1.4 s to 0.6 s for `"heat"` and 1.4 s
#' to 0.6 s for `"rwr"`.
#'
#' Two things are worth knowing about where the remaining time goes. The kernel
#' is built once per call, and for `"laplacian"` that factorization is itself
#' expensive -- 2.8 s on the graph above -- so at small `n_permutations` most of
#' the runtime is a fixed cost that no batching touches. The saving is
#' proportional to `n_permutations`, so it is substantial at the default of
#' 1000 and negligible at 50. If you intend to run several seed sets against one
#' graph, build the kernel once with [prepare_diffusion()] instead.
#'
#' Batching is numerically free for `"laplacian"` and `"heat"`, which produce
#' bit-identical results either way. `"rwr"` is the exception: the iteration
#' stops when the slowest column of a batch has converged, so the other columns
#' are carried slightly past their own stopping point. They end up marginally
#' *closer* to the fixed point, and the difference is bounded by the convergence
#' tolerance, but `"rwr"` scores are therefore reproducible to that tolerance
#' rather than to the bit.
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
#' `null = "degree_matched"` removes that confound by grouping candidates into
#' adaptive degree bins and replacing each real seed with a vertex drawn
#' uniformly from its own bin. The resulting p-values answer "are these seeds
#' unusually well *placed*?" rather than "are these seeds unusually well
#' *connected*?", which is almost always the intended question.
#'
#' How closely the match can be made is limited by the graph. If the seed set is
#' extreme -- the six highest-degree vertices of a scale-free network, say --
#' there are simply no other vertices of comparable degree to permute with, and
#' no permutation null can fully correct for degree. The function warns when the
#' closest available candidates are still far from the seeds, so the residual
#' confound is visible rather than assumed away.
#'
#' **Matching more tightly is not better.** There is a real trade-off here, and
#' it runs the opposite way to intuition. Narrowing the bins improves the degree
#' match but shrinks the set of vertices the null ever places a seed on, and
#' each permutation draws `length(seed_nodes)` of them. Once the bins approach
#' the size of the seed set, consecutive permutations are mostly the same
#' vertices: a node far from all of them scores near zero in every permutation
#' and is then reported as significant whatever its real score. The null becomes
#' *anti*-conservative -- which defeats the purpose, since the degree-matched
#' null exists to be the more conservative of the two.
#'
#' The default therefore scales the pool with the seed count rather than fixing
#' it: `10 * length(seed_nodes)`. That multiplier was calibrated on five
#' synthetic scale-free graphs rather than derived, and should be treated as a
#' working default rather than a settled one.
#'
#' One consequence is worth stating, because it is easy to mistake for a bug.
#' The degree-matched null is **not** uniformly more conservative than the
#' uniform null, node count for node count. Because its bins restrict which
#' vertices are ever seeded, it draws each of them more often -- on a
#' 150-vertex graph with 5 seeds, up to 10 of 60 permutations against 5 for the
#' uniform null. Each vertex is tested only against the permutations that did
#' not seed it, so that narrower support means a larger share of a vertex's own
#' null is excluded, and what gets excluded is its upper tail. The remaining
#' null is smaller and the p-value lower. Both nulls stay valid; they are
#' answering slightly different conditional questions, and the matched one is
#' conditioned more strongly. Narrowing `match_pool` strengthens that
#' conditioning, which is a second reason not to narrow it.
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
                                           match_pool = NULL) {

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

  # 1. Build the diffusion kernel, once.
  #
  # Built here, before the real scores, rather than after them. The real-score
  # call used to run first and without `precompute`, so network_diffusion() built
  # a kernel of its own and threw it away, and the identical kernel was then
  # built again for the permutations. For "laplacian" that is a second sparse
  # Cholesky factorization of the whole graph: 2.8 s of the 6.8 s this function
  # took on a 6000-vertex, 78000-edge graph went on a factorization whose result
  # was discarded. The cost is per call rather than per permutation, so unlike
  # the batching below it is invisible in any benchmark that scales
  # n_permutations up -- and it dominates at small ones.
  precomp <- prepare_diffusion(graph = graph,
                               method = method,
                               alpha = alpha, t = t, restart_prob = restart_prob,
                               normalize = normalize,
                               weights = weights,
                               weight_type = weight_type)

  # 2. Compute real diffusion scores
  real_scores_df <- network_diffusion(graph, seed_nodes, method = method,
                                      alpha = alpha, t = t,
                                      restart_prob = restart_prob,
                                      normalize = normalize,
                                      precompute = precomp,
                                      weights = weights,
                                      weight_type = weight_type,
                                      seed_weights = seed_weights)
  # Pinned to graph order, not to score order. network_diffusion() sorts its
  # output by descending score, so the real scores and every permutation come
  # back in a *different* vertex order. Everything downstream is matched by
  # position, so the alignment has to happen here and in the permutation lambda,
  # once each, before any of it is combined.
  real_scores <- setNames(real_scores_df$score, real_scores_df$node)[all_nodes]

  if (verbose) {
    # class(plan()) is e.g. c("FutureStrategy", "sequential", "uniprocess", ...);
    # element 2 is the strategy name ("sequential", "multisession", ...).
    plan_name <- class(future::plan())[2]
    message(sprintf("Running %d permutations with %d random seed nodes each (future plan: %s)...",
                    n_permutations, n_seeds, plan_name))
  }

  # 3. Run permutations under whatever future plan the caller has set.
  #
  # This deliberately does NOT call future::plan(). A package must not change the
  # user's plan: doing so overrides their choice globally, and the "multisession"
  # workers it starts are never shut down, leaving socket connections open. That
  # surfaces as "checking examples ... ERROR / connections left open" under
  # R CMD check. The default plan is sequential; see @details for opting in to
  # parallel execution.

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

  # Draw every permuted seed set up front, one future_lapply element per
  # permutation with `future.seed = seed`. The batching below changes how the
  # diffusions are *solved*, not which seeds are drawn: keeping the draw in this
  # shape means each permutation still consumes its own L'Ecuyer stream, so the
  # permuted sets -- and therefore the p-values -- are identical to those the
  # per-permutation implementation produced for the same `seed`.
  perm_seed_sets <- future.apply::future_lapply(
    seq_len(n_permutations), function(i) sampler(), future.seed = seed
  )

  # Solve the permutations in batches, as a multi-RHS problem rather than one
  # diffusion per permutation. A sparse Cholesky solve has a fixed overhead per
  # call that does not grow with the number of right-hand sides, so 200 separate
  # solves pay it 200 times; the heat expansion and the random walk amortize the
  # sparse traversal over the columns instead of repeating it per vector.
  #
  # Measured on a 6000-vertex, 78000-edge graph with 200 permutations, the
  # permutation solves alone went from 0.75 s to 0.26 s for "laplacian", 1.01 s
  # to 0.46 s for "heat" and 0.70 s to 0.60 s for "rwr". Note what that implies
  # about the *function*: at 200 permutations most of its runtime is the
  # one-time Cholesky factorization in prepare_diffusion() (2.8 s on this
  # graph), which no amount of batching touches. The batching win is
  # proportional to n_permutations, so it matters at the documented default of
  # 1000 and more at 5000, and is close to invisible at 50.
  #
  # The batches are then distributed with future_lapply, so a caller who has set
  # a plan still gets parallelism -- across batches now rather than across
  # permutations. No RNG is consumed here, the seed sets having already been
  # drawn, hence `future.seed = FALSE`.
  n_nodes <- length(all_nodes)
  batch_cols <- diffusion_batch_size(n_nodes, n_permutations)
  batches <- split(seq_len(n_permutations),
                   ceiling(seq_len(n_permutations) / batch_cols))

  perm_results <- future.apply::future_lapply(batches, function(idx) {
    # Rows in graph order by construction. The per-permutation implementation
    # had to re-align each result by name, because network_diffusion() sorts its
    # output by descending score and so returned every permutation in a
    # different vertex order; binding those by position gave each node the
    # distribution of the k-th *ranked* null score. Building the seed matrix
    # directly removes that step rather than repeating it -- there is no sort to
    # undo -- but the row order is still asserted below, since a matrix whose
    # rows have been permuted looks exactly as well-formed as one whose rows
    # have not.
    F0 <- matrix(
      vapply(perm_seed_sets[idx],
             function(ps) make_seed_vector(all_nodes, ps, seed_weights),
             numeric(n_nodes)),
      nrow = n_nodes
    )
    as.matrix(diffuse_columns(F0, precomp, method = method, t = t,
                              restart_prob = restart_prob))
  }, future.seed = FALSE)

  # 4. One column per permutation, every column in `all_nodes` order, as is
  # `real_scores`, so binding by position is correct.
  null_scores_mat <- do.call(cbind, perm_results)
  rownames(null_scores_mat) <- all_nodes

  # 5. Which permutations seeded which vertex.
  #
  # A permutation in which vertex v was itself a permuted seed contributes an
  # enormous score at v -- for the same reason the real seeds top the real run --
  # and that permutation therefore says nothing about whether v is unusually well
  # placed relative to the real seeds. Leaving it in the null inflates the upper
  # tail of v's own null distribution and makes v's p-value conservative,
  # systematically so for the low-degree vertices that get drawn most often.
  # Each vertex is now tested only against the permutations that did not seed it.
  #
  # Built as a logical matrix in `all_nodes` x permutation order, matching
  # null_scores_mat exactly, so the mask and the scores cannot drift apart.
  seeded <- matrix(FALSE, nrow = n_nodes, ncol = n_permutations,
                   dimnames = list(all_nodes, NULL))
  for (i in seq_len(n_permutations)) {
    seeded[perm_seed_sets[[i]], i] <- TRUE
  }

  # 6. Per-node null summaries and p-values.
  #
  # row-wise rather than via split(): split(m, row(m)) allocates an integer
  # matrix the size of the data and then one vector per vertex, which on a
  # 6000-vertex graph is 6000 allocations to compute 6000 numbers.
  stats_rows <- lapply(seq_len(n_nodes), function(i) {
    keep <- !seeded[i, ]
    nulls <- null_scores_mat[i, keep]
    nulls <- nulls[is.finite(nulls)]
    n_valid <- length(nulls)

    if (n_valid == 0L) {
      # Every permutation seeded this vertex, so it has no null to be tested
      # against. Reported as NA rather than as a p-value computed from nothing.
      return(c(p = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
               z = NA_real_, n_valid = 0))
    }

    mu <- mean(nulls)
    sdev <- stats::sd(nulls)
    c(
      # One-sided, upper tail: the question is whether the real score is higher
      # than the null, never lower. The +1/+1 convention keeps p away from zero,
      # so the smallest attainable value is 1 / (n_valid + 1).
      p = (sum(nulls >= real_scores[[i]]) + 1) / (n_valid + 1),
      null_mean = mu,
      null_sd = sdev,
      z = if (is.finite(sdev) && sdev > 0) {
        (real_scores[[i]] - mu) / sdev
      } else {
        NA_real_
      },
      n_valid = n_valid
    )
  })
  stats_mat <- do.call(rbind, stats_rows)

  # 7. Return result
  result <- tibble::tibble(
    node = names(real_scores),
    score = as.numeric(real_scores),
    p_empirical = stats_mat[, "p"],
    null_mean = stats_mat[, "null_mean"],
    null_sd = stats_mat[, "null_sd"],
    z = stats_mat[, "z"],
    n_valid = as.integer(stats_mat[, "n_valid"])
  )

  # Ordered by p, then by descending z. The p-value floor is 1 / (n_valid + 1),
  # so with the default 1000 permutations every vertex above the threshold piles
  # up at 0.000999 and the ordering within that pile would otherwise be whatever
  # order the vertices happen to sit in. The z-score breaks the tie by effect
  # size, which is the only thing left that distinguishes them.
  result <- result[order(result$p_empirical, -result$z, na.last = TRUE), ]
  return(result)
}

#' Columns per batch for the permutation null
#'
#' Internal. Chooses how many permuted seed vectors to solve in one multi-RHS
#' call. Two constraints pull in opposite directions, so neither a fixed number
#' nor "all of them at once" is right:
#'
#' * The intermediate is dense and `n x batch`, so the batch has to be bounded
#'   by the vertex count rather than fixed. A 6000-vertex graph solved 1000
#'   permutations at once would allocate 48 MB per intermediate, and the heat
#'   expansion holds three of them at a time.
#' * Parallelism now happens across batches. One batch means one worker does
#'   everything however many the caller started, so the batch count is capped to
#'   give every worker something to do.
#'
#' @param n Vertex count.
#' @param n_permutations Total permutations to solve.
#' @param budget Target size of one dense intermediate, in matrix entries.
#'   `1e6` doubles is about 8 MB.
#'
#' @return A single positive integer, never larger than `n_permutations`.
#'
#' @keywords internal
#' @noRd
diffusion_batch_size <- function(n, n_permutations, budget = 1e6) {

  cols <- max(1L, as.integer(floor(budget / max(n, 1L))))
  cols <- min(cols, n_permutations)

  workers <- tryCatch(future::nbrOfWorkers(), error = function(e) 1L)
  if (length(workers) == 1L && is.finite(workers) && workers > 1L) {
    cols <- min(cols, as.integer(ceiling(n_permutations / workers)))
  }

  max(1L, as.integer(cols))
}

#' Build a degree-matched seed sampler
#'
#' Internal. Returns a closure that draws a permuted seed set whose degrees match
#' the real seed set's, by replacing each real seed with a randomly chosen
#' candidate of similar degree.
#'
#' Adaptive degree bins, not per-seed nearest-degree pools. The per-seed version
#' this replaces matched degree tightly and was, as a null, close to useless: it
#' gave each seed the `match_pool` nearest candidates by `order()`, which breaks
#' ties by vertex index, so seeds of similar degree received near-identical
#' pools and then drew from them without replacement. Across 200 permutations of
#' a 5-seed set on a 150-vertex scale-free graph the sampler touched **10 of 145
#' candidates**. A null assembled from ten vertices has almost no variance, so
#' nearly any real score exceeded nearly all of it, and the p-values were
#' anti-conservative -- the degree-matched null reported three times as many
#' significant nodes as the uniform one it is supposed to be *more* conservative
#' than. That defect was invisible until the permutation rows were aligned
#' correctly, because the p-values were scrambled anyway.
#'
#' Binning restores the variance. Both naive schemes really are too loose, for
#' the reasons the previous implementation recorded -- quantiles of degree value
#' collapse on a heavy tail, quantiles of degree rank put a 1 and a 17 in one
#' bin -- so neither is used. Instead the bins are grown adaptively: walk the
#' candidates in degree order and close a bin only once it holds `match_pool`
#' members *and* the degree is about to change. Bins are therefore as narrow as
#' the degree distribution allows while still being wide enough to sample from,
#' which on a scale-free graph means the degree-1 vertices form one large bin and
#' the hubs form bins only as wide as they have to be.
#'
#' Three further properties matter, each of which the previous version got
#' wrong: ties are broken at random rather than by vertex index; seeds are
#' processed in a random order each draw, so the first seed is not always the one
#' with the widest remaining choice; and the drawn set is returned in the
#' *original* seed order, because `network_diffusion_with_pvalues()` reuses
#' `seed_weights` against it positionally.
#'
#' @param graph An `igraph` object.
#' @param seed_nodes The real seed names.
#' @param candidates Vertex names eligible as permuted seeds.
#' @param match_pool The minimum number of candidates in each degree bin, and so
#'   the smallest pool any seed may be drawn from. Small values match degree
#'   tightly but give the null little variance; `Inf` puts every candidate in one
#'   bin and so reduces exactly to the uniform null.
#'
#' @return A function of no arguments returning a character vector of seed names.
#'
#' @keywords internal
#' @noRd
make_degree_matched_sampler <- function(graph, seed_nodes, candidates,
                                        match_pool) {

  n_seeds <- length(seed_nodes)

  if (!is.null(match_pool) &&
      (!is.numeric(match_pool) || length(match_pool) != 1 || match_pool < 1)) {
    stop("'match_pool' must be a single number of at least 1.", call. = FALSE)
  }

  # Each permutation consumes n_seeds members of a bin, so a bin barely larger
  # than the seed set yields permutations that are mostly the same vertices.
  # The null is then narrow in *placement* as well as in degree: a node far
  # from every bin member scores zero in every permutation and so comes out
  # significant whatever its real score. Measured on five scale-free graphs,
  # a floor of 10x the seed count is where the degree-matched null becomes
  # reliably more conservative than the uniform one -- which is what it is for.
  # Wider than about 20x and the bins span too much degree to be a match at all.
  POOL_PER_SEED <- 10L
  auto_pool <- min(max(POOL_PER_SEED * n_seeds, 10L), length(candidates))

  if (is.null(match_pool)) {
    match_pool <- auto_pool
  } else if (match_pool < POOL_PER_SEED * n_seeds &&
             match_pool < length(candidates)) {
    warning(sprintf(
      paste0("match_pool = %s is small for %d seeds: each permutation would ",
             "draw %d of as few as %s candidates, leaving the null too narrow ",
             "and the p-values anti-conservative. Consider match_pool = %d or ",
             "more, or match_pool = NULL to choose it automatically."),
      format(match_pool), n_seeds, n_seeds, format(match_pool),
      auto_pool
    ), call. = FALSE)
  }

  deg <- igraph::degree(graph)
  names(deg) <- igraph::V(graph)$name

  # When every candidate is eligible for every seed, this *is* the uniform null.
  # Returning that sampler rather than an equivalent-in-distribution one makes
  # the equivalence exact, which is what the documentation claims, and avoids
  # warning below about a match that was never requested.
  if (match_pool >= length(candidates)) {
    return(function() sample(candidates, n_seeds))
  }

  cand_deg <- deg[candidates]
  min_bin <- max(1L, as.integer(min(match_pool, length(candidates))))

  # --- Adaptive degree bins over the candidates ---
  #
  # Walk the candidates in degree order, closing a bin once it holds min_bin
  # members and the next degree differs. Splitting a run of equal degrees would
  # be an arbitrary cut, and a bin below min_bin is too small to be sampled
  # from, so both constraints have to hold before a bin is closed.
  ord <- order(cand_deg)
  sorted_names <- candidates[ord]
  sorted_deg <- unname(cand_deg[ord])

  bin_id <- integer(length(sorted_deg))
  current <- 1L
  filled <- 0L
  for (i in seq_along(sorted_deg)) {
    bin_id[i] <- current
    filled <- filled + 1L
    if (filled >= min_bin && i < length(sorted_deg) &&
        sorted_deg[i + 1L] != sorted_deg[i]) {
      current <- current + 1L
      filled <- 0L
    }
  }
  # A short trailing bin is merged back rather than left underfilled, so the
  # occupancy floor holds for every bin and not merely for most of them.
  if (current > 1L && filled < min_bin) {
    bin_id[bin_id == current] <- current - 1L
  }

  members <- split(sorted_names, factor(bin_id, levels = sort(unique(bin_id))))
  bin_lo <- vapply(members, function(m) min(deg[m]), numeric(1))
  bin_hi <- vapply(members, function(m) max(deg[m]), numeric(1))
  bin_mean <- vapply(members, function(m) mean(deg[m]), numeric(1))

  # Each seed takes the bin whose degree range contains it, or the nearest bin
  # when its degree lies outside every range -- which is what happens to a hub
  # more connected than any available candidate.
  seed_bin <- vapply(seed_nodes, function(s) {
    d <- deg[[s]]
    gap <- pmax(bin_lo - d, 0) + pmax(d - bin_hi, 0)
    which.min(gap)
  }, integer(1))

  # How well the match can possibly do is limited by what the graph contains: an
  # extreme seed set has no comparable vertices to be replaced by. Say so once,
  # rather than leaving the caller to assume the correction was exact.
  #
  # Judged against the closest candidates that *exist*, not against the bin
  # actually used. The bins are deliberately wider than the tightest possible
  # match, so measuring the bin would report that deliberate widening as though
  # it were a property of the graph -- and would fire on nearly every call.
  seed_deg <- deg[seed_nodes]
  n_closest <- min(5L, length(cand_deg))
  best <- vapply(seed_nodes, function(s) {
    mean(sort(cand_deg[order(abs(cand_deg - deg[[s]]))][seq_len(n_closest)]))
  }, numeric(1))

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
    out <- character(n_seeds)

    # Random processing order. With a fixed one the first seed always picks from
    # a full bin and the last from whatever is left, which biases the match by
    # position in the seed vector.
    for (i in sample.int(n_seeds)) {
      pool <- setdiff(members[[seed_bin[[i]]]], drawn)

      if (length(pool) == 0) {
        # The bin is exhausted because other seeds took its members. Widen to
        # every remaining candidate and keep the closest in degree, so the match
        # degrades gracefully instead of failing.
        others <- setdiff(candidates, drawn)
        if (length(others) == 0) {
          stop("Not enough candidate vertices to build a degree-matched ",
               "permutation. Use null = \"uniform\", or reduce the seed set.",
               call. = FALSE)
        }
        gaps <- abs(deg[others] - deg[[seed_nodes[[i]]]])
        pool <- others[gaps == min(gaps)]
      }

      # sample.int() over the pool, so equally good candidates are equally
      # likely. The previous version resolved ties through order(), which is
      # stable and therefore always returned the lowest vertex index.
      pick <- pool[[sample.int(length(pool), 1)]]
      drawn <- c(drawn, pick)
      out[[i]] <- pick
    }

    # Original seed order, not draw order: network_diffusion_with_pvalues()
    # reuses `seed_weights` against this vector positionally, so permuting it
    # would attach each weight to the wrong replacement.
    out
  }
}
