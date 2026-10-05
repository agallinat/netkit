#' Compute a Table of Node-Level Centrality Metrics
#'
#' Computes several per-node centrality and position metrics in one call and
#' returns them as a tibble, with every metric also attached to the graph as a
#' vertex attribute.
#'
#' This is the node-level counterpart to [summarize_graph_metrics()], which
#' describes the graph as a whole. Having both means a ranking of nodes no longer
#' has to be stitched together from individual \pkg{igraph} calls, and because the
#' returned graph carries the metrics as vertex attributes, the result feeds
#' directly into [plot_Net()] (`color = "pagerank"`) and into
#' [robustness_analysis()] (`removal_strategy = "pagerank"`), which already
#' accepts any numeric vertex attribute as a removal priority.
#'
#' @param graph An `igraph` object representing the network to analyze, or a data
#'   frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param metrics Character vector of metrics to compute. Any subset of
#'   `"degree"`, `"strength"`, `"betweenness"`, `"closeness"`, `"harmonic"`,
#'   `"eigenvector"`, `"pagerank"`, `"coreness"`, `"clustering"`, `"constraint"`
#'   and `"eccentricity"`. Defaults to all but `"closeness"` -- see Details.
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#' @param normalized Logical. If `TRUE` (default), metrics that \pkg{igraph} can
#'   normalize to `[0, 1]` are normalized, which makes them comparable across
#'   graphs of different size.
#' @param mode For directed graphs, whether to count `"all"` (default), `"in"` or
#'   `"out"` edges. Ignored for undirected graphs.
#' @param max_nodes Integer. Above this vertex count, the metrics whose cost is
#'   O(VE) -- `"betweenness"`, `"closeness"`, `"harmonic"` and `"eccentricity"` --
#'   are skipped with a warning and returned as `NA`. Set to `Inf` to compute them
#'   regardless. Default is `5000`, matching the threshold
#'   [summarize_graph_metrics()] already uses for betweenness.
#' @param plot Logical. Whether to build a diagnostic plot. Default is `TRUE`.
#' @param plot_type Either `"correlation"` (default), a Spearman correlation
#'   heatmap of the computed metrics, or `"ranking"`, a faceted bar chart of the
#'   `top_n` highest-scoring nodes per metric.
#' @param top_n Integer. Number of nodes shown per facet when
#'   `plot_type = "ranking"`. Default is `15`.
#' @param label.size Numeric. Base font size for plot text. Default is `12`.
#'
#' @return A list with:
#' \describe{
#'   \item{`plot`}{A `ggplot2` object, or `NULL` when `plot = FALSE`. The element
#'     is always present, so the return shape does not depend on the arguments.}
#'   \item{`result`}{A tibble with one row per vertex: `node` followed by one
#'     column per requested metric, in the order given by `metrics`.}
#'   \item{`graph`}{The input graph with each metric attached as a vertex
#'     attribute of the same name.}
#'   \item{`method`}{A human-readable description of what was computed.}
#' }
#'
#' @details
#' The metrics are:
#'
#' \describe{
#'   \item{`degree`}{Number of incident edges. Always the unweighted count.}
#'   \item{`strength`}{Sum of incident edge weights; equal to `degree` when
#'     unweighted.}
#'   \item{`betweenness`}{Fraction of shortest paths passing through the node.
#'     Uses edge *costs* when weighted.}
#'   \item{`closeness`}{Reciprocal of the mean distance to all *reachable*
#'     nodes. Not comparable across components -- see below.}
#'   \item{`harmonic`}{Sum of the reciprocal distances. Unlike `closeness` this
#'     is comparable across components, because an unreachable pair contributes
#'     \eqn{1/\infty = 0} rather than being excluded from the average.}
#'   \item{`eigenvector`}{Leading eigenvector of the adjacency matrix: a node is
#'     central if its neighbors are.}
#'   \item{`pagerank`}{Stationary distribution of a random surfer. Sums to 1.}
#'   \item{`coreness`}{Largest *k* for which the node belongs to the k-core.}
#'   \item{`clustering`}{Local transitivity: how interconnected the node's
#'     neighborhood is. `NaN` for nodes of degree below 2, which have no
#'     neighbor pairs.}
#'   \item{`constraint`}{Burt's constraint: how much the node's connections are
#'     concentrated within a single group. Low constraint marks a broker.}
#'   \item{`eccentricity`}{Distance to the furthest reachable node.}
#' }
#'
#' `"closeness"` is omitted from the defaults because it is not comparable across
#' components, and most real networks are disconnected. \pkg{igraph} averages the
#' distance over reachable vertices only, so a node in a *small* isolated
#' component scores *higher* than an equally well-placed node in the giant
#' component -- everything in its component is nearby. Ranking by closeness on a
#' fragmented graph therefore puts the periphery on top. `"harmonic"` answers the
#' same question without that failure mode, because unreachable pairs contribute
#' zero rather than being dropped from the denominator. Request `"closeness"`
#' explicitly if the graph is connected or you want the per-component reading.
#'
#' The default plot is a correlation heatmap rather than a ranking because the
#' most common mistake with a table like this is to treat the metrics as
#' independent evidence. On many networks betweenness, eigenvector centrality and
#' PageRank correlate with degree above 0.9, so a node that looks important by
#' four measures may be important by one. The heatmap makes that visible before
#' the ranking is interpreted.
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @references
#' Freeman, L. C. (1978). Centrality in social networks: conceptual
#' clarification. *Social Networks*, 1(3), 215-239.
#' \doi{10.1016/0378-8733(78)90021-7}
#'
#' Burt, R. S. (1992). *Structural Holes: The Social Structure of Competition*.
#' Harvard University Press.
#'
#' Marchiori, M., & Latora, V. (2000). Harmony in the small-world.
#' *Physica A*, 285(3-4), 539-546. \doi{10.1016/S0378-4371(00)00311-3}
#'
#' @seealso [summarize_graph_metrics()] for the graph-level counterpart,
#'   [find_hubs()] for thresholded classification.
#'
#' @examples
#' g <- igraph::sample_pa(60, power = 1.5, directed = FALSE)
#' igraph::V(g)$name <- paste0("n", seq_len(igraph::vcount(g)))
#'
#' res <- node_metrics(g, plot = FALSE)
#' head(res$result)
#'
#' # Every metric is on the returned graph, so it chains straight onward.
#' igraph::vertex_attr_names(res$graph)
#' rob <- robustness_analysis(res$graph, removal_strategy = "pagerank",
#'                            steps = 10, plot = FALSE)
#' rob$auc$lcc_size
#'
#' # The default plot shows how far the metrics actually disagree.
#' p <- node_metrics(g, metrics = c("degree", "betweenness", "pagerank"))$plot
#'
#' @importFrom igraph degree strength betweenness closeness harmonic_centrality
#' @importFrom igraph eigen_centrality page_rank coreness transitivity constraint
#' @importFrom igraph eccentricity vcount vertex_attr vertex_attr<- is_directed
#' @importFrom tibble tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 labs
#' @importFrom ggplot2 theme_minimal theme element_blank geom_col facet_wrap
#' @importFrom ggplot2 coord_flip scale_x_discrete element_text
#' @importFrom stats cor reorder
#'
#' @export
node_metrics <- function(graph,
                         metrics = c("degree", "strength", "betweenness",
                                     "harmonic", "eigenvector", "pagerank",
                                     "coreness", "clustering", "constraint",
                                     "eccentricity"),
                         weights = NULL,
                         weight_type = c("strength", "distance"),
                         normalized = TRUE,
                         mode = c("all", "in", "out"),
                         max_nodes = 5000,
                         plot = TRUE,
                         plot_type = c("correlation", "ranking"),
                         top_n = 15,
                         label.size = 12) {

  all_metrics <- c("degree", "strength", "betweenness", "closeness", "harmonic",
                   "eigenvector", "pagerank", "coreness", "clustering",
                   "constraint", "eccentricity")

  metrics <- match.arg(metrics, choices = all_metrics, several.ok = TRUE)
  mode <- match.arg(mode)
  plot_type <- match.arg(plot_type)

  # Results are keyed by vertex name, so names must exist.
  graph <- as_netkit_graph(graph, backfill_names = TRUE)
  w <- as_netkit_weights(graph, weights, weight_type)

  n <- vcount(graph)
  node_names <- vertex_attr(graph, "name")

  # Metrics that need all-pairs shortest paths are O(VE) and will hang on a large
  # graph. Skipping them loudly, and returning NA rather than nothing, keeps the
  # column set stable regardless of graph size -- the same reason
  # summarize_graph_metrics() reports NA for Avg_betweenness above this bound.
  costly <- c("betweenness", "closeness", "harmonic", "eccentricity")
  skip <- character(0)
  if (n > max_nodes) {
    skip <- intersect(metrics, costly)
    if (length(skip) > 0) {
      warning(sprintf(
        paste0("Graph has %d vertices (max_nodes = %s), so %s %s skipped and ",
               "returned as NA. Raise max_nodes to compute %s anyway."),
        n, format(max_nodes), paste(skip, collapse = ", "),
        if (length(skip) == 1) "was" else "were",
        if (length(skip) == 1) "it" else "them"
      ), call. = FALSE)
    }
  }

  compute <- function(metric) {
    if (metric %in% skip) {
      return(rep(NA_real_, n))
    }
    switch(
      metric,
      degree = as.numeric(igraph::degree(graph, mode = mode,
                                         normalized = normalized)),
      # strength has no `normalized` argument; it is a sum, not a count, so
      # there is no natural maximum to divide by.
      strength = if (w$weighted) {
        as.numeric(igraph::strength(graph, mode = mode, weights = w$strength))
      } else {
        as.numeric(igraph::degree(graph, mode = mode))
      },
      betweenness = as.numeric(igraph::betweenness(graph, weights = w$distance,
                                                   normalized = normalized)),
      closeness = as.numeric(igraph::closeness(graph, mode = mode,
                                               weights = w$distance,
                                               normalized = normalized)),
      harmonic = as.numeric(igraph::harmonic_centrality(graph, mode = mode,
                                                        weights = w$distance,
                                                        normalized = normalized)),
      eigenvector = as.numeric(igraph::eigen_centrality(
        graph, directed = is_directed(graph), weights = w$strength
      )$vector),
      pagerank = as.numeric(igraph::page_rank(graph, weights = w$strength)$vector),
      coreness = as.numeric(igraph::coreness(graph, mode = mode)),
      # type = "local" returns NaN for degree < 2, which is correct: such a node
      # has no pairs of neighbors and therefore no neighborhood to be clustered.
      clustering = as.numeric(igraph::transitivity(graph, type = "local")),
      constraint = as.numeric(igraph::constraint(graph, weights = w$strength)),
      eccentricity = as.numeric(igraph::eccentricity(graph, mode = mode,
                                                     weights = w$distance))
    )
  }

  values <- lapply(metrics, compute)
  names(values) <- metrics

  result <- tibble::as_tibble(c(list(node = node_names), values))

  # Attach to the graph so the metrics chain onward.
  for (metric in metrics) {
    igraph::vertex_attr(graph, metric) <- values[[metric]]
  }

  p <- NULL
  if (plot) {
    p <- if (plot_type == "correlation") {
      node_metrics_corr_plot(result, metrics, label.size)
    } else {
      node_metrics_ranking_plot(result, metrics, top_n, label.size)
    }
  }

  list(
    plot = p,
    result = result,
    graph = graph,
    method = paste0(
      "Node metrics (", paste(metrics, collapse = ", "), "); ",
      describe_weights(w), "; ",
      if (normalized) "normalized" else "unnormalized",
      if (is_directed(graph)) paste0("; mode = ", mode) else "",
      if (length(skip) > 0) {
        paste0("; skipped above max_nodes: ", paste(skip, collapse = ", "))
      } else {
        ""
      }
    )
  )
}

#' Spearman correlation heatmap of node metrics
#'
#' Internal helper for `node_metrics()`.
#'
#' @param result The metric tibble.
#' @param metrics Metric names, in display order.
#' @param label.size Base font size.
#'
#' @return A `ggplot` object, or `NULL` when fewer than two metrics vary.
#'
#' @keywords internal
#' @noRd
node_metrics_corr_plot <- function(result, metrics, label.size) {

  mat <- as.matrix(result[, metrics, drop = FALSE])

  # A constant column has zero variance, so its correlation is undefined and
  # cor() would warn. Drop it rather than emit a warning the caller cannot act on.
  varies <- vapply(seq_along(metrics), function(j) {
    v <- mat[, j]
    sum(is.finite(v)) > 1 && stats::sd(v[is.finite(v)]) > 0
  }, logical(1))

  keep <- metrics[varies]
  if (length(keep) < 2) {
    return(NULL)
  }

  # Spearman rather than Pearson: these metrics are heavy-tailed and what matters
  # is whether they rank the nodes the same way, not whether they are linearly
  # related. "pairwise.complete.obs" handles the NaNs that local clustering
  # legitimately produces for low-degree nodes.
  cm <- stats::cor(mat[, keep, drop = FALSE], method = "spearman",
                   use = "pairwise.complete.obs")

  df <- data.frame(
    metric_x = factor(rep(keep, times = length(keep)), levels = keep),
    metric_y = factor(rep(keep, each = length(keep)), levels = rev(keep)),
    correlation = as.numeric(cm)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = metric_x, y = metric_y,
                                   fill = correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", correlation)),
                       size = label.size / 4) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white",
                                 high = "#b2182b", midpoint = 0,
                                 limits = c(-1, 1)) +
    ggplot2::labs(
      x = NULL, y = NULL, fill = "Spearman",
      title = "Rank correlation between node metrics",
      subtitle = "Highly correlated metrics are not independent evidence"
    ) +
    ggplot2::theme_minimal(base_size = label.size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Faceted top-n ranking of node metrics
#'
#' Internal helper for `node_metrics()`.
#'
#' @param result The metric tibble.
#' @param metrics Metric names, in display order.
#' @param top_n Nodes per facet.
#' @param label.size Base font size.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
node_metrics_ranking_plot <- function(result, metrics, top_n, label.size) {

  pieces <- lapply(metrics, function(metric) {
    v <- result[[metric]]
    ok <- is.finite(v)
    if (!any(ok)) {
      return(NULL)
    }
    ord <- order(v[ok], decreasing = TRUE)
    idx <- which(ok)[utils::head(ord, top_n)]
    data.frame(
      metric = metric,
      node = result$node[idx],
      value = v[idx],
      # A single ordering cannot serve every facet, since each ranks different
      # nodes. Rank within the facet and let facet_wrap use a free scale.
      rank = seq_along(idx),
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, pieces[!vapply(pieces, is.null, logical(1))])
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  # The same node can appear in several facets at different ranks, so one global
  # factor ordering cannot serve them all. Qualify each bar by its facet to make
  # the levels unique, order those levels by (metric, rank), then strip the
  # qualifier back off at draw time so the axis still reads as node names.
  key <- paste(df$node, df$metric, sep = "@@")
  df$label <- factor(key, levels = rev(key[order(df$metric, df$rank)]))

  ggplot2::ggplot(df, ggplot2::aes(x = label, y = value)) +
    ggplot2::geom_col(fill = "#2166ac", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ metric, scales = "free") +
    ggplot2::scale_x_discrete(labels = function(x) sub("^.*", "", x)) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0("Top ", top_n, " nodes per metric")) +
    ggplot2::theme_minimal(base_size = label.size)
}
