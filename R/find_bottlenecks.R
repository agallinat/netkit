#' Identify and Bottleneck Nodes in a Network
#'
#' Identifies bottleneck nodes in an \code{igraph} network as those with low degree
#' and high betweenness centrality. The function supports both standardized (z-score)
#' and quantile-based thresholding. Optionally, it produces a 2D scatter plot with
#' bottlenecks highlighted.
#'
#' @param graph An \code{igraph} object representing the network to analyze or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param method Character. Method to define bottlenecks: \code{"absolute"}
#'   (default), \code{"zscore"} or \code{"quantile"}. See "Choosing a method"
#'   below -- the default changed in 0.1.2.
#' @param max_degree Numeric. Upper bound on the number of neighbors a bottleneck
#'   may have (only used if \code{method = "absolute"}). Default is \code{3}.
#' @param degree_threshold Numeric. Upper threshold for standardized degree (only used if \code{method = "zscore"}).
#' @param betweenness_threshold Numeric. Lower threshold for standardized betweenness (used in both methods).
#' @param degree_quantile Numeric between 0 and 1. Quantile threshold for degree (used if \code{method = "quantile"}).
#' @param betweenness_quantile Numeric between 0 and 1. Quantile threshold for betweenness (used if \code{method = "quantile"}).
#' @param log_transform Logical. If \code{TRUE}, applies \code{log1p} transformation to degree and betweenness.
#' @param plot Logical. If \code{TRUE}, generates a plot of degree vs. betweenness highlighting bottlenecks.
#' @param focus_color Character. Color to display in the focus area of the plot (bottlenecks region).
#' @param bottleneck_names Logical. If \code{TRUE}, labels bottleneck nodes on the plot.
#' @param bottleneck_cex Numeric. Font size scaling for bottleneck labels on the plot.
#' @param gg_extra List. Additional user-defined layers for the returned ggplot.
#' eg. list(ylim(-2,2), theme_bw(), theme(legend.position = "none"))
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{plot}}{A scatter plot of degree vs. betweenness highlighting
#'     bottlenecks, or \code{NULL} when \code{plot = FALSE}. The element is always
#'     present, so the return shape does not depend on the arguments. Note this is a
#'     \code{ggExtraPlot} (a \pkg{grid} gtable), not a \code{ggplot}: the
#'     scatterplot is passed through \code{\link[ggExtra]{ggMarginal}()} to add the
#'     marginal histograms, which assembles it, so it can be printed but not
#'     extended with \code{+}. To customize it, pass layers through
#'     \code{gg_extra}; those are added before the assembly step.}
#'   \item{\code{method}}{A message describing the method and thresholds used.}
#'   \item{\code{result}}{A \code{tibble} with node name, degree, betweenness, transformed metrics, and bottleneck status.}
#'   \item{\code{graph}}{The original graph with a new vertex attribute \code{is_bottleneck} (logical).}
#' }
#'
#' @section Empty selections:
#'
#' A bottleneck must satisfy two criteria at once -- low degree *and* high
#' betweenness -- and the conjunction is frequently empty for reasons that have
#' nothing to do with the thresholds. In a dense graph there are low-degree
#' vertices and there are high-betweenness vertices, but rarely one that is
#' both: with many alternative routes available, a poorly connected vertex lies
#' on few shortest paths. That is a real answer to the question, and the
#' function now says so rather than returning an empty table in silence. The
#' warning distinguishes it from the other two causes -- a threshold no vertex
#' can reach, and a metric with no variance at all.
#'
#' @section Choosing a method:
#'
#' The default is `"absolute"`, which was `"zscore"` through 0.1.1. The two axes
#' of a bottleneck are not the same kind of statement, and treating them as
#' though they were is what made the old default unreachable:
#'
#' * **Low degree is absolute.** A bottleneck is poorly connected -- it has few
#'   neighbors, full stop. `"zscore"` instead asked for low degree to be a
#'   *statistical outlier*, with `degree_threshold = -1` requiring a degree a
#'   standard deviation below the mean. On any right-skewed degree distribution
#'   the left tail is short: degree-1 vertices are the single most common kind of
#'   vertex in a scale-free graph, so they are the *mode*, not an anomaly, and
#'   their z-score is around -0.6. The default therefore selected nothing at
#'   every size from n = 60 to n = 2000 -- not because those graphs have no
#'   bottlenecks, but because "unusually poorly connected" is close to a
#'   contradiction in a graph where being poorly connected is the norm.
#' * **High betweenness is relative.** Carrying a lot of flow only means
#'   something by comparison with the rest of the graph, so this axis stays a
#'   quantile (`betweenness_quantile`, default the 75th percentile).
#'
#' So `"absolute"` selects vertices with at most `max_degree` neighbors whose
#' betweenness is in the top `1 - betweenness_quantile` of the vertices that lie
#' on any shortest path at all. Two details matter:
#'
#' * The degree axis is the neighbor **count**, even when `weights` is supplied,
#'   and is neither log-transformed nor standardized. An absolute threshold needs
#'   a unit that means the same thing in every graph, and only the neighbor count
#'   qualifies: a strength cap of 3 means nothing without knowing what the
#'   weights measure. Betweenness still honours `weights`, so connectivity is
#'   counted while flow is weighted.
#' * The betweenness quantile is computed over the vertices with **nonzero**
#'   betweenness. Over the full distribution it would be zero on any sparse
#'   graph, for the reason given below.
#'
#' `max_degree = 3` is a default, not a law. It is low enough to mean "poorly
#' connected" in most networks and high enough to admit a genuine multi-edged
#' bridge; raise it on a dense graph, where no vertex may be that sparse, and
#' note that a vertex of degree 1 can never be selected whatever it is set to,
#' since a leaf lies on no shortest path between other vertices.
#'
#' `"zscore"` and `"quantile"` are unchanged and still available. See the
#' "Choosing a method" section of [find_hubs()] for the trade-off between those
#' two, which applies identically here -- with the caveat that the `"zscore"`
#' degree default is the one described above and will usually select nothing.
#'
#' A vertex with betweenness of exactly zero is never selected, whatever the
#' thresholds. A bottleneck lies on shortest paths, and such a vertex lies on
#' none. This is a necessary guard rather than a refinement: betweenness is
#' heavily zero-inflated on a sparse graph -- 87% of the vertices of a
#' 2000-vertex scale-free graph have betweenness 0 -- which puts the 75th
#' percentile *at* zero, and an inclusive `>=` comparison against zero then
#' admits every one of them. The default quantile call used to return 1731 of
#' 2000 vertices as bottlenecks, every one of them with betweenness exactly
#' zero. No choice of quantile defends against that, because the quantile of a
#' mostly-zero distribution is zero.
#'
#' When a metric takes the same value at every vertex the standardization is
#' undefined. `degree_metric` and `betweenness_metric` are then `NA`, and
#' `is_bottleneck` is `FALSE`: no vertex is unusual when none differs from any
#' other.
#'
#' @examples
#' # Two dense clusters joined by a two-vertex bridge. The bridge vertices have
#' # low degree and carry every path between the clusters, which is exactly the
#' # low-degree, high-betweenness signature a bottleneck is defined by -- a
#' # scale-free graph has no such vertex, so the defaults find nothing there.
#' g <- igraph::disjoint_union(igraph::make_full_graph(6),
#'                             igraph::make_full_graph(6))
#' igraph::V(g)$name <- c(paste0("L", 1:6), paste0("R", 1:6))
#' g <- igraph::add_vertices(g, 2, name = c("b1", "b2"))
#' g <- igraph::add_edges(g, c("L1", "b1", "b1", "b2", "b2", "R1"))
#'
#' res <- find_bottlenecks(g, plot = FALSE)
#' res$method
#' res$result[res$result$is_bottleneck, c("node", "degree", "betweenness")]
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
find_bottlenecks <- function(graph,
                             method = c("absolute", "zscore", "quantile"),
                             max_degree = 3,
                             degree_threshold = -1,
                             betweenness_threshold = 1,
                             degree_quantile = 0.25,
                             betweenness_quantile = 0.75,
                             log_transform = TRUE,
                             plot = TRUE,
                             focus_color = "skyblue",
                             bottleneck_names = TRUE,
                             bottleneck_cex = 3,
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
  deg_count <- degree(graph)
  deg <- if (w$weighted) igraph::strength(graph, weights = w$strength) else deg_count
  btw <- betweenness(graph, normalized = TRUE, weights = w$distance)

  # Optional log transform
  deg_val <- if (log_transform) log1p(deg) else deg
  btw_val <- if (log_transform) log1p(btw) else btw

  deg_degenerate <- FALSE
  btw_degenerate <- FALSE

  if (method == "absolute") {
    if (!is.numeric(max_degree) || length(max_degree) != 1L ||
        !is.finite(max_degree) || max_degree < 1) {
      stop("`max_degree` must be a single finite number of at least 1.")
    }

    # The degree axis is the neighbor *count*, even under weights, and is
    # deliberately not log-transformed or standardized. An absolute threshold
    # needs a unit that means the same thing in every graph, and the number of
    # neighbors is the only such unit available: a strength cap of 3 says nothing
    # without knowing what the weights measure, and log1p(degree) <= 3 is not a
    # statement anyone intends. Betweenness still honours `weights`, so
    # connectivity is counted while flow is weighted.
    deg_score <- deg_count
    btw_score <- btw_val

    # The betweenness cutoff is taken over the vertices that lie on at least one
    # shortest path, not over all of them. On a sparse graph most vertices have
    # betweenness exactly zero, so a quantile of the full distribution is itself
    # zero -- the trap documented below, which no choice of quantile escapes.
    # Restricted to the positive values, the cutoff is a real quantile of the
    # distribution of actual flow.
    positive <- btw > 0
    btw_pos <- btw_score[positive]
    btw_cutoff <- if (length(btw_pos)) {
      quantile(btw_pos, betweenness_quantile, na.rm = TRUE)
    } else {
      Inf
    }

    # "High betweenness" is a comparison, so it is undefined when the vertices
    # carrying flow all carry the same amount. Without this guard a quantile
    # cutoff admits all of them -- every vertex of a ring is at the 75th
    # percentile of a constant -- which is precisely the criticism of taking a
    # fixed fraction of the graph regardless of its structure. One flow-carrying
    # vertex is not degenerate: it is unambiguously the maximum.
    s_pos <- if (length(btw_pos) > 1L) stats::sd(btw_pos) else NA_real_
    btw_degenerate <- length(btw_pos) > 1L && (!is.finite(s_pos) || s_pos == 0)

    degree_threshold <- max_degree
    betweenness_threshold <- btw_cutoff
    is_bottleneck <- if (btw_degenerate) {
      rep(FALSE, length(deg_score))
    } else {
      deg_score <= max_degree & btw_score >= btw_cutoff
    }
  } else if (method == "zscore") {
    # Via standardize_node_metric() rather than scale() directly: on a graph
    # where every vertex has the same degree, scale() returns NaN throughout and
    # the comparisons below are then NA, not FALSE.
    deg_std <- standardize_node_metric(deg_val)
    btw_std <- standardize_node_metric(btw_val)
    deg_score <- deg_std$z
    btw_score <- btw_std$z
    deg_degenerate <- deg_std$degenerate
    btw_degenerate <- btw_std$degenerate
    is_bottleneck <- deg_score <= degree_threshold & btw_score >= betweenness_threshold
  } else if (method == "quantile") {
    deg_score <- deg_val
    btw_score <- btw_val
    deg_cutoff <- quantile(deg_score, degree_quantile, na.rm = TRUE)
    btw_cutoff <- quantile(btw_score, betweenness_quantile, na.rm = TRUE)
    degree_threshold <- deg_cutoff
    betweenness_threshold <- btw_cutoff
    is_bottleneck <- deg_score <= degree_threshold & btw_score >= betweenness_threshold
  }

  # A vertex whose metric is undefined is not a bottleneck. Leaving the NA in
  # place wrote an all-NA vertex attribute onto the graph, which every
  # downstream filter then silently dropped.
  is_bottleneck[is.na(is_bottleneck)] <- FALSE

  # A bottleneck lies on shortest paths. A vertex with betweenness of exactly
  # zero lies on none, so it cannot be one whatever the thresholds say -- and
  # under `method = "quantile"` they frequently did say so. Betweenness is
  # heavily zero-inflated on a sparse graph: 87% of the vertices of a
  # 2000-vertex scale-free graph have betweenness 0, which puts the 75th
  # percentile *at* zero, and `>=` then admits every one of them. The default
  # call returned 1731 of 2000 vertices as bottlenecks, every single one with
  # betweenness exactly 0. A quantile cutoff cannot defend against this on its
  # own: the problem is that the quantile of a mostly-zero distribution is zero.
  is_bottleneck <- is_bottleneck & btw > 0

  # An empty selection used to be returned in silence. Here the conjunction is
  # the usual culprit rather than the thresholds: a dense graph has low-degree
  # vertices and high-betweenness vertices but rarely a vertex that is both, so
  # the warning distinguishes that from an unreachable threshold.
  if (!any(is_bottleneck)) {
    criterion <- if (method == "zscore") {
      c("degree_threshold", "betweenness_threshold")
    } else if (method == "absolute") {
      c("max_degree",
        sprintf("the betweenness %s quantile of the vertices with nonzero betweenness",
                format(betweenness_quantile)))
    } else {
      c(sprintf("the degree %s quantile", format(degree_quantile)),
        sprintf("the betweenness %s quantile", format(betweenness_quantile)))
    }
    warn_empty_selection(
      "find_bottlenecks", "bottleneck nodes", method,
      list(
        list(metric = criterion[1], label = "degree metric", direction = "<=",
             threshold = degree_threshold,
             attained = suppressWarnings(min(deg_score, na.rm = TRUE)),
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
    is_bottleneck = is_bottleneck
  )

  # Create vertex attribute 'bottleneck'
  vertex_attr(graph, "is_bottleneck") <- is_bottleneck

  # Plot if requested
  if (plot) {
    df <- result
    p <- ggplot(df, aes(x = degree_metric, y = betweenness_metric)) +
      annotate("rect",
               xmin = -Inf, xmax = degree_threshold,
               ymin = betweenness_threshold, ymax = Inf,
               alpha = 0.2, fill = focus_color) +
      geom_point(aes(color = is_bottleneck), size = 2, alpha = 0.8) +
      scale_color_manual(values = c("gray80", "#e41a1c")) +
      geom_vline(xintercept = degree_threshold, linetype = "dashed", color = focus_color) +
      geom_hline(yintercept = betweenness_threshold, linetype = "dashed", color = focus_color) +
      geom_text_repel(data = filter(df, is_bottleneck),
                      aes(label = node), color = "#e41a1c",
                      point.padding =  3,
                      size = ifelse(bottleneck_names, bottleneck_cex, 0)) +  # Tamaño de los nombres.
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = "Degree Metric", y = "Betweenness Metric", color = "Bottleneck")

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
        method = paste("Bottlenecks identified by method:", method,
                       "with Degree metric threshold =", degree_threshold,
                       "and Betweenness metric threshold =", betweenness_threshold,
                     paste0("(", describe_weights(w), ")")),
        result = result,
        graph = graph
      ), "find_bottlenecks"))

}
