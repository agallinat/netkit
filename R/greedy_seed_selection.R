#' Greedy Seed Node Selection to Maximize Diffusion Toward Target Nodes
#'
#' This function implements a greedy algorithm to select a set of `k` seed nodes from a candidate list
#' such that the resulting diffusion signal on a specified set of target nodes is maximized.
#' It supports multiple diffusion models (Laplacian, Heat Kernel, Random Walk with Restart).
#'
#' @param graph An \code{igraph} object representing the network to analyze or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param target_nodes A character vector of node names to prioritize for receiving the diffusion signal.
#' @param candidate_nodes Optional character vector of eligible nodes to consider as seeds.
#' If NULL (default), all non-target nodes are used.
#' @param method Diffusion method to use. One of `"laplacian"`, `"heat"`, or `"rwr"`.
#' @param alpha Scaling parameter for the Laplacian diffusion method (default = 0.7).
#' @param t Time parameter for the Heat Kernel diffusion (default = 1).
#' @param restart_prob Restart probability for the Random Walk with Restart (default = 0.3).
#' @param k Number of seed nodes to select (default = 5).
#' @param normalize Logical; whether to symmetrically normalize the adjacency
#'   matrix before diffusing, passed through to [network_diffusion()]
#'   (default = TRUE). This controls the diffusion kernel, *not* the target
#'   score: it was previously documented as normalizing the cumulative score
#'   over the target nodes, which it has never done.
#' @param plot Logical; whether to plot the progression of target diffusion score as seeds are selected (default = TRUE).
#'
#' @param objective What the greedy search maximizes:
#'   \describe{
#'     \item{`"total"`}{The default, and the historical behaviour: the summed
#'       diffusion score over `target_nodes`. Be aware that this makes the
#'       greedy search redundant -- see "Why the greedy search is not greedy".}
#'     \item{`"coverage"`}{A saturating objective, which makes the search
#'       genuinely greedy: signal delivered to a target that is already well
#'       served counts for progressively less, so the selection spreads across
#'       the target set instead of concentrating on whichever target is easiest
#'       to reach.}
#'   }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{selected_seeds}{Character vector of selected seed node names, in the
#'     order they were chosen.}
#'   \item{final_target_score}{Final cumulative diffusion score over the target nodes.}
#'   \item{scores_at_each_step}{Vector of target scores at each greedy selection
#'     step. Always the summed diffusion score over the targets, whichever
#'     `objective` was used, so the two objectives' selections are comparable on
#'     one scale.}
#'   \item{objective_at_each_step}{The value of the objective actually being
#'     maximized at each step. Identical to `scores_at_each_step` under
#'     `objective = "total"`; under `"coverage"` it is an effective number of
#'     targets reached, bounded by `length(target_nodes)`.}
#'   \item{result}{A tibble with one row per selected seed: `step`, `node`,
#'     `degree`, `gain` (the marginal gain that pick bought), `gain_z` (that
#'     gain standardized against candidates of comparable degree -- see "Reading
#'     `gain_z`"), `cumulative_score` and `objective`.}
#'   \item{plot}{A `ggplot2` object showing the progression, or `NULL`.}
#'   \item{method}{A description of the objective and diffusion used, including
#'     whether the selection reduces to a top-k ranking.}
#' }
#'
#' @details
#' This function uses a stepwise greedy heuristic. At each step it evaluates the
#' marginal gain in the objective from adding each candidate to the current seed
#' set and selects the one with the highest gain, repeated for `k` steps.
#'
#' The target score is the sum of the diffusion scores over `target_nodes`. It
#' is a total, never an average, and `normalize` does not change that -- it
#' changes the kernel the scores come from.
#'
#' @section Why the greedy search is not greedy:
#'
#' Under `objective = "total"` this function cannot do better than ranking the
#' candidates once, and it is worth being explicit about why.
#'
#' Diffusion is a *linear* operator. The signal a seed set delivers to a target
#' is therefore exactly the sum of what each member delivers on its own, so the
#' total target score is a modular set function: the gain from adding a
#' candidate does not depend on which candidates were already chosen. Every step
#' of the search consequently ranks the remaining candidates in the same order
#' as the first step did, and the result is precisely the `k` candidates with
#' the largest single-seed contributions. A greedy search is only worth running
#' against an objective with *interacting* terms, and this one has none.
#'
#' `objective = "coverage"` supplies the interaction. Each target's accumulated
#' signal is passed through a saturating transform before being summed, so the
#' second seed to serve an already-well-covered target adds little and the
#' search is pushed toward targets nothing has reached yet. That function is
#' submodular rather than modular, which is the condition under which greedy
#' selection is both non-trivial and known to be near-optimal.
#'
#' @section What gets selected, and why it may not be what you want:
#'
#' With `normalize = TRUE` the kernel is symmetrically normalized, so an edge
#' into a vertex of degree \eqn{d} carries weight \eqn{1/\sqrt{d}}. A vertex of
#' degree 1 adjacent to a target is therefore the most efficient single place to
#' inject signal that exists in the graph: it has nowhere else to send it.
#'
#' The consequence is that this function will often return leaves hanging off
#' the targets. That is the correct answer to the question as posed and usually
#' not the answer that was wanted, since a degree-1 vertex in a biological
#' network is more often a gap in the annotation than a controllable
#' intervention point. The function warns when more than half the selection has
#' degree 1. To get something actionable, restrict `candidate_nodes` to vertices
#' you could actually perturb, or set `normalize = FALSE` to stop rewarding low
#' degree.
#'
#' @section Reading `gain_z`:
#'
#' Because the raw `gain` is largely a restatement of degree, `result` also
#' reports `gain_z`: the pick's gain standardized against the gains of the
#' candidates closest to it in degree. It separates the two reasons a vertex can
#' look good.
#'
#' A `gain_z` near zero means the vertex is typical of its degree class -- it
#' was selected *because of* its degree, and any other vertex of the same degree
#' in the same position would have served as well. A large `gain_z` means the
#' vertex outperforms others as poorly (or as well) connected as itself, which
#' is a statement about where it sits in the graph rather than about how many
#' edges it has. When a selection is full of leaves, `gain_z` is what tells you
#' whether any of them is individually worth attention.
#'
#' Two details to read it correctly. The comparison pool is the candidates still
#' *available* at that step, not the original candidate set, so `gain_z` answers
#' "how unusual was this pick among the alternatives at the moment it was made".
#' A consequence is that it can climb across steps while `gain` stays flat: once
#' the other strong leaves have been taken, the leaf picked fourth really does
#' stand out further among what remains. That is not increasing importance.
#'
#' And it is `NA` only when every comparable candidate has an identical gain,
#' which needs a vertex-transitive arrangement such as a complete graph.
#' Equal *degrees* are not enough and are in fact the best case: on a ring every
#' candidate has degree 2, so the pool is the whole candidate set and `gain_z`
#' measures position with no degree confound left in it at all.
#'
#' @section Performance:
#'
#' Every diffusion is computed once, in a single batched pass over
#' `candidate_nodes`, and the search is then arithmetic on the resulting
#' `length(target_nodes)` x `length(candidate_nodes)` matrix. The previous
#' implementation ran a full diffusion for every candidate at every step --
#' `k * length(candidate_nodes)` diffusions, recomputing at step 3 what it had
#' already computed at step 1.
#'
#' On a 2000-vertex graph with 4 targets and `k = 8`, that is 15968 diffusions
#' and 15.7 s, against 0.37 s now, for the same selection and the same scores to
#' within 3e-16.
#'
#' One consequence is worth noting for `method = "rwr"`, whose iterative solve
#' is only additive up to its convergence tolerance: the selected *set* is
#' unchanged, but the order within it may differ from older versions by
#' differences of around 1e-9, which were arbitrary to begin with.
#'
#' @seealso \code{\link{network_diffusion}}, \code{\link{network_diffusion_with_pvalues}}
#'
#' @examples
#' g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
#' igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
#' target <- c("1", "2", "3")
#'
#' res <- greedy_seed_selection(g, target_nodes = target, k = 3, plot = FALSE)
#' res$selected_seeds
#' res$scores_at_each_step
#'
#' # The default objective is additive, so that selection is just the top 3
#' # single-seed contributions. `"coverage"` saturates each target's share, which
#' # spreads the selection over the target set instead.
#' cov <- greedy_seed_selection(g, target_nodes = target, k = 3,
#'                              objective = "coverage", plot = FALSE)
#' cov$selected_seeds
#' cov$objective_at_each_step   # bounded by length(target_nodes)
#'
#' @importFrom igraph is_igraph vertex_attr vertex_attr<- vcount
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @importFrom progress progress_bar
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous labs theme_minimal
#'
#' @export
greedy_seed_selection <- function(graph,
                                  target_nodes,
                                  candidate_nodes = NULL,
                                  method = c("laplacian", "heat", "rwr"),
                                  alpha = 0.7, t = 1, restart_prob = 0.3,
                                  k = 5,
                                  normalize = TRUE,
                                  plot = TRUE,
                                  weights = NULL,
                                  weight_type = c("strength", "distance"),
                                  objective = c("total", "coverage")) {

  # --- Input validation ---
  method <- match.arg(method)
  weight_type <- match.arg(weight_type)
  objective <- match.arg(objective)

  graph <- as_netkit_graph(graph, backfill_names = TRUE)

  all_nodes <- igraph::V(graph)$name
  target_nodes <- intersect(target_nodes, all_nodes)

  if (length(target_nodes) == 0) {
    stop("No target nodes found in the graph.")
  }

  if (is.null(candidate_nodes)) {
    candidate_nodes <- setdiff(all_nodes, target_nodes)
  }

  if (length(candidate_nodes) == 0) {
    stop("No candidate nodes available to select from.", call. = FALSE)
  }

  k <- min(k, length(candidate_nodes))

  precomp <- prepare_diffusion(
    graph = graph,
    method = method,
    alpha = alpha,
    t = t,
    restart_prob = restart_prob,
    normalize = normalize,
    weights = weights,
    weight_type = weight_type
  )

  # --- Single-seed contributions to each target ---
  #
  # Every diffusion this function needs is computed here, once. Diffusion is a
  # linear operator, so the signal a seed set delivers to a target is the sum of
  # what its members deliver individually; the whole search therefore reduces to
  # arithmetic on one |targets| x |candidates| matrix. The previous
  # implementation ran a full diffusion per candidate per step --
  # `k * length(candidate_nodes)` of them, recomputing at step 3 exactly what it
  # had already computed at step 1.
  contrib <- target_contributions(graph, target_nodes, candidate_nodes, precomp,
                                  method = method, t = t,
                                  restart_prob = restart_prob)

  # --- Greedy selection ---
  #
  # Run as an actual greedy loop for both objectives rather than sorting for the
  # linear one. The two coincide under `"total"` -- that is the point of the
  # warning below -- but making the code rely on that would bake the equivalence
  # in as an assumption instead of leaving it a property the tests can check.
  cand_degree <- igraph::degree(graph, v = candidate_nodes)
  sel <- greedy_select(contrib, k = k, objective = objective,
                       degree = cand_degree)

  selected <- candidate_nodes[sel$order]
  best_scores <- sel$total
  objective_scores <- sel$objective

  # --- Degree of the chosen seeds ---
  #
  # Under symmetric normalization an edge into a degree-1 vertex carries weight
  # 1/sqrt(d_target), the largest any edge to that target can carry, so a leaf
  # hanging off a target is the single most efficient place to inject signal.
  # That is arithmetically correct and usually useless: such a vertex is
  # typically an artefact of incomplete annotation rather than a controllable
  # intervention point. Worth saying out loud, because nothing in the returned
  # scores hints at it.
  sel_degree <- igraph::degree(graph, v = selected)
  if (mean(sel_degree == 1) > 0.5) {
    warning(sprintf(
      paste0("%d of %d selected seeds have degree 1. Diffusion normalization ",
             "makes a leaf adjacent to a target the most efficient single ",
             "injection point, so the selection is dominated by low-degree ",
             "satellites rather than by well-connected regulators. Restrict ",
             "`candidate_nodes` to the vertices you could actually ",
             "intervene on, or set normalize = FALSE, if that is not what you ",
             "want."),
      sum(sel_degree == 1), length(selected)
    ), call. = FALSE)
  }

  # --- Optional plotting ---
  if (plot && length(best_scores) > 1) {
    df_plot <- data.frame(step = seq_along(best_scores), target_score = best_scores)

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = step, y = target_score)) +
      ggplot2::geom_line(color = "#e41a1c") +
      ggplot2::geom_point(color = "#e41a1c") +
      ggplot2::scale_x_continuous(breaks = 1:k) +
      ggplot2::labs(
        title = "Target Scores by Greedy Seed Selection Step",
        x = "Step",
        y = ifelse(normalize, "Normalized Target Score", "Target Score")
      ) +
      ggplot2::theme_minimal(base_size = 14)

  } else {
    p <- NULL
  }

  # --- Return results ---
  result <- tibble::tibble(
    step = seq_along(selected),
    node = selected,
    degree = as.integer(unname(cand_degree[sel$order])),
    gain = sel$gain,
    gain_z = sel$gain_z,
    cumulative_score = best_scores,
    objective = objective_scores
  )

  return(as_netkit_result(list(
    selected_seeds = selected,
    final_target_score = utils::tail(best_scores, 1),
    scores_at_each_step = best_scores,
    objective_at_each_step = objective_scores,
    result = result,
    plot = p,
    method = paste0(
      "Greedy seed selection maximizing the '", objective, "' objective over ",
      length(target_nodes), " target nodes (", method, " diffusion",
      if (normalize) ", normalized" else "", "); ",
      if (objective == "total") {
        paste0("the total is additive, so this selection is identical to the ",
               "top-", length(selected),
               " single-seed contributions -- see the 'Why the greedy search ",
               "is not greedy' section")
      } else {
        "saturating per-target gains, so the selection spreads across targets"
      }
    )
  ), "greedy_seed_selection"))
}

#' Single-seed contribution of every candidate to every target
#'
#' Internal. Returns the `length(target_nodes)` x `length(candidate_nodes)`
#' matrix whose `[v, c]` entry is the diffusion score at target `v` when `c` is
#' the only seed.
#'
#' Computed in batches of candidate columns through `diffuse_columns()`, so the
#' whole matrix costs roughly one diffusion pass over the candidate set rather
#' than one diffusion per candidate. Only the target rows are retained, so the
#' memory held is `|targets| x |candidates|` however large the graph is -- the
#' `n x batch` intermediate is transient.
#'
#' @param graph An `igraph` object with named vertices.
#' @param target_nodes,candidate_nodes Character vectors of vertex names.
#' @param precomp A kernel from [prepare_diffusion()].
#' @param method,t,restart_prob Diffusion parameters, as in [network_diffusion()].
#'
#' @return A numeric matrix with one row per target and one column per
#'   candidate, in the order the names were given.
#'
#' @keywords internal
#' @noRd
target_contributions <- function(graph, target_nodes, candidate_nodes, precomp,
                                 method, t = 1, restart_prob = 0.3) {

  all_nodes <- igraph::V(graph)$name
  n <- length(all_nodes)
  target_idx <- match(target_nodes, all_nodes)

  n_cand <- length(candidate_nodes)
  batch <- diffusion_batch_size(n, n_cand)
  batches <- split(seq_len(n_cand), ceiling(seq_len(n_cand) / batch))

  out <- matrix(NA_real_, nrow = length(target_nodes), ncol = n_cand,
                dimnames = list(target_nodes, candidate_nodes))

  # The contribution pass is now the only expensive step, so the progress bar
  # belongs here rather than on the greedy loop it used to track -- that loop is
  # arithmetic on a matrix that already exists and finishes instantly.
  pb <- progress::progress_bar$new(
    format = "  Diffusing candidates [:bar] :percent eta: :eta",
    total = length(batches), clear = FALSE, width = 60
  )

  for (idx in batches) {
    pb$tick()
    F0 <- matrix(
      vapply(candidate_nodes[idx],
             function(cand) make_seed_vector(all_nodes, cand, NULL),
             numeric(n)),
      nrow = n
    )
    scores <- diffuse_columns(F0, precomp, method = method, t = t,
                              restart_prob = restart_prob)
    # drop = FALSE: a single target would otherwise collapse the matrix to a
    # vector and silently transpose the assignment.
    out[, idx] <- as.matrix(scores)[target_idx, , drop = FALSE]
  }

  out
}

#' Greedy maximization over a contribution matrix
#'
#' Internal. Picks `k` columns of `contrib` one at a time, each time taking the
#' column with the largest marginal gain in the chosen objective.
#'
#' The two objectives differ in exactly one respect, and it is the respect that
#' decides whether a greedy search is doing anything:
#'
#' * `"total"` sums the per-target signal. Because diffusion is linear this is a
#'   *modular* set function -- the gain from adding a candidate does not depend
#'   on what is already selected -- so every step ranks the candidates
#'   identically and the result is the top `k` single-seed contributions.
#' * `"coverage"` sums a saturating transform of each target's accumulated
#'   signal, `1 - exp(-x / m_v)`, where `m_v` is the most any single candidate
#'   delivers to target `v`. A concave increasing function of a modular function
#'   is submodular, so gains genuinely diminish: once a target is well served,
#'   further signal to it counts for little and the search moves on to the
#'   targets that are not. The value is bounded by one per target, and so reads
#'   as an effective number of targets reached.
#'
#' `1 - exp(-x/m_v)` rather than `min(x/m_v, 1)`: a hard cap gives exactly zero
#' marginal gain once a target saturates, which leaves the argmax undetermined
#' among candidates that only reach saturated targets. The exponential keeps
#' every gain strictly positive while still diminishing, so the pick is always
#' well defined.
#'
#' Targets that no candidate reaches at all (`m_v == 0`) are dropped from the
#' coverage sum. Keeping them would divide by zero; weighting them would be
#' pretending a disconnected target can be served.
#'
#' @param contrib A targets x candidates matrix of single-seed contributions.
#' @param k Number of columns to select.
#' @param objective Either `"total"` or `"coverage"`.
#' @param degree Candidate degrees, aligned with the columns of `contrib`. Used
#'   only to standardize each pick's gain against candidates of comparable
#'   degree; see `degree_matched_gain_z()`.
#'
#' @return A list with `order` (selected column indices, in selection order),
#'   `total` and `objective` (the summed target score and the objective value
#'   after each step), `gain` (the marginal gain the pick bought) and `gain_z`
#'   (that gain standardized against comparable-degree candidates).
#'
#' @keywords internal
#' @noRd
greedy_select <- function(contrib, k, objective = c("total", "coverage"),
                          degree = NULL) {

  objective <- match.arg(objective)

  n_cand <- ncol(contrib)
  # Per-target reference level. rowMaxs by hand: apply() over a one-row matrix
  # would drop to a vector.
  m <- apply(contrib, 1, max)
  usable <- is.finite(m) & m > 0

  value <- function(acc) {
    if (objective == "total") return(sum(acc))
    if (!any(usable)) return(0)
    sum(1 - exp(-acc[usable] / m[usable]))
  }

  acc <- rep(0, nrow(contrib))
  current <- value(acc)
  remaining <- seq_len(n_cand)
  chosen <- integer(0)
  totals <- numeric(0)
  objectives <- numeric(0)
  gains <- numeric(0)
  gain_zs <- numeric(0)

  for (step in seq_len(k)) {
    step_values <- vapply(remaining, function(j) value(acc + contrib[, j]),
                          numeric(1))
    winner <- which.max(step_values)
    best <- remaining[[winner]]

    # Marginal gains of every candidate still available at this step. The
    # z-score has to be taken against these, not against the single-seed
    # contributions, because under "coverage" the two differ once anything has
    # been selected.
    step_gains <- step_values - current

    gain_zs <- c(gain_zs, if (is.null(degree)) {
      NA_real_
    } else {
      degree_matched_gain_z(step_gains, degree[remaining], winner)
    })
    gains <- c(gains, step_gains[[winner]])

    acc <- acc + contrib[, best]
    current <- value(acc)
    chosen <- c(chosen, best)
    remaining <- setdiff(remaining, best)
    totals <- c(totals, sum(acc))
    objectives <- c(objectives, current)
  }

  list(order = chosen, total = totals, objective = objectives,
       gain = gains, gain_z = gain_zs)
}

#' Standardize a marginal gain against comparable-degree candidates
#'
#' Internal. Returns how unusual a candidate's gain is *for a candidate of its
#' degree*, as a z-score against the gains of the candidates closest to it in
#' degree.
#'
#' This exists because the raw gain is mostly a statement about degree. Under
#' symmetric normalization a leaf adjacent to a target is the most efficient
#' injection point in the graph, so leaves have the largest gains and the
#' selection fills up with them -- but if every leaf has a similar gain, no
#' individual leaf is informative. A gain z-score near zero says "this vertex is
#' typical of its degree class and was picked for its degree"; a large one says
#' "this vertex beats others as poorly connected as itself", which is the
#' claim worth acting on.
#'
#' The comparison pool is every candidate within the smallest degree radius that
#' admits at least `min_pool` of them. Taking a whole radius rather than the
#' `min_pool` nearest avoids resolving ties by vertex index -- the mistake that
#' made the diffusion null's degree-matched sampler degenerate. The candidate
#' itself is included in its own pool, as a z-score within a reference group
#' conventionally is.
#'
#' @param gains Marginal gains of the candidates still available.
#' @param degrees Their degrees, aligned with `gains`.
#' @param i Index into `gains` of the candidate being standardized.
#' @param min_pool Smallest acceptable comparison pool.
#'
#' @return A single numeric, `NA` when the pool has no spread of gains -- which
#'   happens legitimately when every comparable candidate is equivalent.
#'
#' @keywords internal
#' @noRd
degree_matched_gain_z <- function(gains, degrees, i, min_pool = 10L) {

  dist <- abs(degrees - degrees[[i]])
  radii <- sort(unique(dist))
  # Smallest radius admitting min_pool candidates; the largest radius covers
  # everything, so this always resolves.
  counts <- vapply(radii, function(r) sum(dist <= r), integer(1))
  radius <- radii[[which(counts >= min(min_pool, length(gains)))[[1]]]]

  pool <- gains[dist <= radius]
  sdev <- stats::sd(pool)
  if (!is.finite(sdev) || sdev == 0) return(NA_real_)
  (gains[[i]] - mean(pool)) / sdev
}


