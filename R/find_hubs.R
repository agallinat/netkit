#' Identify Hub Nodes in a Network
#'
#' This function identifies hub nodes in an \code{igraph} network based on
#' degree and betweenness centrality using either z-score or quantile thresholds.
#' Optionally, it visualizes the classification using a scatter plot with marginal histograms.
#'
#' @param graph An \code{igraph} object or a data frame containing a symbolic edge list in the
#'   first two columns. Additional columns are considered as edge attributes.
#' @param method Character. Method to identify hubs: \code{"zscore"} (standardized metrics)
#'   or \code{"quantile"} (empirical percentiles).
#' @param degree_threshold Numeric. Threshold for standardized degree (only used if \code{method = "zscore"}).
#' @param betweenness_threshold Numeric. Threshold for standardized betweenness (only used if \code{method = "zscore"}).
#' @param degree_quantile Numeric between 0 and 1. Quantile threshold for degree (used if \code{method = "quantile"}).
#' @param betweenness_quantile Numeric between 0 and 1. Quantile threshold for betweenness (used if \code{method = "quantile"}).
#' @param log_transform Logical. If \code{TRUE}, applies log-transformation to degree and betweenness metrics.
#' @param plot Logical. If \code{TRUE}, generates a plot of the degree vs. betweenness classification.
#' @param focus_color Character. Color to display in the focus area of the plot (hubs region).
#' @param label.size Numeric. Base font size for plot elements. Passed to \code{theme_classic}.
#' @param hub_names Logical. If \code{TRUE}, adds node labels to identified hubs on the plot.
#' @param hub_cex Numeric. Font size scaling factor for hub labels on the plot.
#' @param gg_extra List. Additional user-defined layers for the returned ggplot.
#' eg. list(ylim(-2,2), theme_bw(), theme(legend.position = "none"))
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{plot}}{A scatter plot of degree vs. betweenness with hub nodes
#'     highlighted, or \code{NULL} when \code{plot = FALSE}. The element is always
#'     present, so the return shape does not depend on the arguments. Note this is a
#'     \code{ggExtraPlot} (a \pkg{grid} gtable), not a \code{ggplot}: the
#'     scatterplot is passed through \code{\link[ggExtra]{ggMarginal}()} to add the
#'     marginal histograms, which assembles it, so it can be printed but not
#'     extended with \code{+}. To customize it, pass layers through
#'     \code{gg_extra}; those are added before the assembly step.}
#'   \item{\code{method}}{Description of the method and thresholds used.}
#'   \item{\code{result}}{A \code{tibble} with node name, degree, betweenness, transformed metrics, and hub status.}
#'   \item{\code{graph}}{The original graph with a new vertex attribute \code{is_hub}.}
#' }
#'
#' @section Choosing a method:
#'
#' The two methods ask different questions and fail in opposite directions.
#'
#' `"zscore"` asks whether a vertex is unusual *relative to the spread of this
#' graph*, which is the meaningful definition of a hub, but its thresholds are
#' only reachable if the metric has a tail to reach them. The defaults of 3 and
#' 1 are calibrated for a scale-free degree distribution. On a dense graph -- a
#' thresholded correlation matrix, say -- the degree distribution is
#' near-Gaussian, the largest z-score a few hundred vertices can produce is
#' around 2 to 3, and `degree_threshold = 3` is therefore *unreachable* rather
#' than unmet. The function now warns when that happens instead of returning an
#' empty selection in silence.
#'
#' `"quantile"` asks which vertices are in the top few percent, which is always
#' answerable -- and that is its weakness: it returns a non-empty selection
#' whether or not the graph has any hub structure at all. On a ring, where every
#' vertex is identical, `degree_quantile = 0.95` marks *every* vertex a hub,
#' because every tied value is at or above the 95th percentile. Use it when you
#' want a fixed-size shortlist from a graph you already believe has hubs; use
#' `"zscore"` when the existence of hubs is itself the question.
#'
#' When a metric takes the same value at every vertex the standardization is
#' undefined. `degree_metric` and `betweenness_metric` are then `NA`, which is
#' what an undefined z-score should report, while `is_hub` is `FALSE`: no vertex
#' is unusual when none differs from any other.
#'
#' A vertex with betweenness of exactly zero is never a hub, whatever the degree
#' threshold says. A hub is defined by high degree *and* high betweenness, and
#' zero is not high. On a complete graph no vertex mediates any path, so none is
#' a hub however many edges it carries -- which is the right answer. In practice
#' this guard changes nothing on a scale-free graph, where a high-degree vertex
#' has paths through it; it is there because "in practice" is not a contract.
#' See [find_bottlenecks()], where the same arithmetic had a much larger effect.
#'
#' @examples
#' g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
#' res <- find_hubs(g, method = "quantile", plot = FALSE)
#' res$method
#' head(res$result)
#'
#' # The returned graph carries an `is_hub` vertex attribute, so results chain.
#' table(igraph::V(res$graph)$is_hub)
#'
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. When supplied, degree becomes vertex *strength* and
#'   betweenness is computed on edge *costs*. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @importFrom igraph is_igraph degree betweenness vertex_attr_names vertex_attr vertex_attr<- vcount
#' @importFrom tibble tibble
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes annotate geom_point scale_color_manual geom_vline geom_hline labs theme_classic theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggExtra ggMarginal
#' @importFrom stats quantile
#'
#' @export
find_hubs <- function(graph,
                      method = c("zscore", "quantile"),
                      degree_threshold = 3,
                      betweenness_threshold = 1,
                      degree_quantile = 0.95,
                      betweenness_quantile = 0.95,
                      log_transform = TRUE,
                      plot = TRUE,
                      focus_color = "darkgreen",
                      label.size = 12,
                      hub_names = TRUE,
                      hub_cex = 3,
                      gg_extra = list(),
                      weights = NULL,
                      weight_type = c("strength", "distance")) {

  method <- match.arg(method)

  graph <- as_netkit_graph(graph, backfill_names = TRUE)
  w <- as_netkit_weights(graph, weights, weight_type)

  # Compute degree and betweenness. Under weights these become vertex strength
  # and cost-based betweenness: a "high degree, high betweenness" node in a
  # weighted graph is one with strong attachments on many shortest paths, which
  # is not the same set of nodes as the unweighted counts pick out.
  deg <- if (w$weighted) igraph::strength(graph, weights = w$strength) else degree(graph)
  btw <- betweenness(graph, normalized = TRUE, weights = w$distance)

  # Optional log transform
  deg_val <- if (log_transform) log1p(deg) else deg
  btw_val <- if (log_transform) log1p(btw) else btw

  deg_degenerate <- FALSE
  btw_degenerate <- FALSE

  if (method == "zscore") {
    # Via standardize_node_metric() rather than scale() directly: on a graph
    # where every vertex has the same degree, scale() returns NaN throughout and
    # the comparisons below are then NA, not FALSE.
    deg_std <- standardize_node_metric(deg_val)
    btw_std <- standardize_node_metric(btw_val)
    deg_score <- deg_std$z
    btw_score <- btw_std$z
    deg_degenerate <- deg_std$degenerate
    btw_degenerate <- btw_std$degenerate
    is_hub <- deg_score >= degree_threshold & btw_score >= betweenness_threshold
  } else if (method == "quantile") {
    deg_score <- deg_val
    btw_score <- btw_val
    deg_cutoff <- quantile(deg_score, degree_quantile, na.rm = TRUE)
    btw_cutoff <- quantile(btw_score, betweenness_quantile, na.rm = TRUE)
    degree_threshold <- deg_cutoff
    betweenness_threshold <- btw_cutoff
    is_hub <- deg_score >= degree_threshold & btw_score >= betweenness_threshold
  }

  # A vertex whose metric is undefined is not a hub. Leaving the NA in place
  # wrote an all-NA vertex attribute onto the graph, which every downstream
  # filter then silently dropped.
  is_hub[is.na(is_hub)] <- FALSE

  # A hub is defined by high degree *and* high betweenness, and zero is not
  # high. Betweenness is heavily zero-inflated on a sparse graph, so a quantile
  # cutoff can land exactly on zero and then admit every vertex on that axis,
  # leaving the selection resting on degree alone while still claiming both.
  # See find_bottlenecks(), where the same arithmetic returned 1731 of 2000
  # vertices. It bites far less often here, because a high-degree vertex usually
  # has paths through it -- but "usually" is not a contract.
  is_hub <- is_hub & btw > 0

  # An empty selection used to be returned in silence, which on a dense
  # correlation graph is the normal outcome of the z-score defaults: the degree
  # distribution has no tail to reach a z of 3, so the answer is "the threshold
  # is unreachable", not "this graph has no hubs".
  if (!any(is_hub)) {
    criterion <- if (method == "zscore") {
      c("degree_threshold", "betweenness_threshold")
    } else {
      c(sprintf("the degree %s quantile", format(degree_quantile)),
        sprintf("the betweenness %s quantile", format(betweenness_quantile)))
    }
    warn_empty_selection(
      "find_hubs", "hub nodes", method,
      list(
        list(metric = criterion[1], label = "degree metric", direction = ">=",
             threshold = degree_threshold,
             attained = suppressWarnings(max(deg_score, na.rm = TRUE)),
             degenerate = deg_degenerate),
        list(metric = criterion[2], label = "betweenness metric",
             direction = ">=", threshold = betweenness_threshold,
             attained = suppressWarnings(max(btw_score, na.rm = TRUE)),
             degenerate = btw_degenerate)
      )
    )
  }

  # Compile result
  # `degree` always holds the unweighted count and `strength` the weighted sum,
  # so neither column changes meaning with the arguments -- `degree_metric` is the
  # one that reflects whichever was actually thresholded. Overloading `degree` to
  # hold strengths would have made the table silently incomparable between calls.
  result <- tibble(
    node = vertex_attr(graph, "name"),
    degree = igraph::degree(graph),
    strength = if (w$weighted) igraph::strength(graph, weights = w$strength) else igraph::degree(graph),
    betweenness = btw,
    degree_metric = deg_score,
    betweenness_metric = btw_score,
    is_hub = is_hub
  )

  # Create vertex attribute 'bottleneck'
  vertex_attr(graph, "is_hub") <- is_hub

  # Plot if requested
  if (plot) {
    df <- result
    p <- ggplot(df, aes(x = degree_metric, y = betweenness_metric)) +
      annotate("rect",
               xmin = degree_threshold, xmax = Inf,
               ymin = betweenness_threshold, ymax = Inf,
               alpha = 0.2, fill = focus_color) +
      geom_point(aes(color = is_hub), size = 2, alpha = 0.8) +
      scale_color_manual(values = c("gray80", "#e41a1c")) +
      geom_vline(xintercept = degree_threshold, linetype = "dashed", color = focus_color) +
      geom_hline(yintercept = betweenness_threshold, linetype = "dashed", color = focus_color) +
      geom_text_repel(data = filter(df, is_hub),
                      aes(label = node), color = "#e41a1c",
                      point.padding =  3,
                      size = ifelse(hub_names, hub_cex, 0)) +  # Names size,
      theme_classic(base_size = label.size) +
      theme(legend.position = "bottom") +
      labs(x = "Degree Metric", y = "Betweenness Metric", color = "Hub nodes")

    # ---- add user‑supplied ggplot layers/themes/etc ----
    if (length(gg_extra)) {
      for (layer in gg_extra) p <- p + layer
    }

    p <- ggExtra::ggMarginal(p, type = "histogram", bins = 30, fill = "lightgray")

  } else {

    p <- NULL

  }

  return(
    as_netkit_result(list(
      plot = p,
      method = paste("Hub nodes identified by method:", method,
                     "with Degree metric threshold =", degree_threshold,
                     "and Betweenness metric threshold =", betweenness_threshold,
                     paste0("(", describe_weights(w), ")")),
      result = result,
      graph = graph),
      "find_hubs")
    )

}
