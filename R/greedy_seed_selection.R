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
#' @param normalize Logical; whether to normalize the cumulative diffusion score over the target nodes (default = TRUE).
#' @param plot Logical; whether to plot the progression of target diffusion score as seeds are selected (default = TRUE).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{selected}{Character vector of selected seed node names.}
#'   \item{final_target_score}{Final cumulative (or normalized) diffusion score over the target nodes.}
#'   \item{scores_at_each_step}{Vector of target scores at each greedy selection step.}
#' }
#'
#' @details
#' This function uses a stepwise greedy heuristic. At each step, it evaluates the marginal gain in target score from
#' adding each candidate node to the current seed set, and selects the one with the highest gain. This is repeated for `k` steps.
#'
#' Note: The score of each candidate at every step is recomputed via a full diffusion run, making the function computationally
#' intensive for large graphs or large `k`.
#'
#' The target score can be interpreted as either the total or average diffusion signal received by the target nodes.
#' Normalization helps scale results across networks of different sizes.
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
                                  weight_type = c("strength", "distance")) {

  # --- Input validation ---
  method <- match.arg(method)
  weight_type <- match.arg(weight_type)

  graph <- as_netkit_graph(graph, backfill_names = TRUE)

  all_nodes <- igraph::V(graph)$name
  target_nodes <- intersect(target_nodes, all_nodes)

  if (length(target_nodes) == 0) {
    stop("No target nodes found in the graph.")
  }

  if (is.null(candidate_nodes)) {
    candidate_nodes <- setdiff(all_nodes, target_nodes)
  }

  selected <- character(0)
  best_scores <- numeric(0)

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

  # --- Progress bar ---
  pb <- progress::progress_bar$new(
    format = "  Selecting seeds [:bar] :percent eta: :eta",
    total = k, clear = FALSE, width = 60
  )

  # --- Helper function for computing target score ---
  compute_target_score <- function(seed_set) {
    scores <- network_diffusion(
      graph = graph,
      seed_nodes = seed_set,
      method = method,
      alpha = alpha,
      t = t,
      restart_prob = restart_prob,
      normalize = normalize,
      precompute = precomp,
      weights = weights,
      weight_type = weight_type
    )
    sum(scores$score[scores$node %in% target_nodes])
  }

  # --- Greedy selection loop ---
  for (i in seq_len(k)) {
    scores_list <- lapply(candidate_nodes, function(node) {
        compute_target_score(c(selected, node))
      })

    best_node <- candidate_nodes[[which.max(scores_list)]]
    selected <- c(selected, best_node)
    candidate_nodes <- setdiff(candidate_nodes, best_node)
    best_scores <- c(best_scores, max(unlist(scores_list)))

    pb$tick()
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
  return(as_netkit_result(list(
    selected_seeds = selected,
    final_target_score = utils::tail(best_scores, 1),
    scores_at_each_step = best_scores,
    plot = p
  ), "greedy_seed_selection"))
}


