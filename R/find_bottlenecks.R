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
#' @param method Character. Method to define bottlenecks: \code{"zscore"} or \code{"quantile"}.
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
#'   \item{\code{method}}{A message describing the method and thresholds used.}
#'   \item{\code{result}}{A \code{tibble} with node name, degree, betweenness, transformed metrics, and bottleneck status.}
#'   \item{\code{graph}}{The original graph with a new vertex attribute \code{bottleneck} (logical).}
#' }
#' If \code{plot = TRUE}, a scatter plot of degree vs. betweenness is displayed, highlighting bottlenecks.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- sample_pa(100)
#' find_bottlenecks(g, method = "quantile", plot = TRUE)
#' }
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
                             method = c("zscore", "quantile"),
                             degree_threshold = -1,
                             betweenness_threshold = 1,
                             degree_quantile = 0.25,
                             betweenness_quantile = 0.75,
                             log_transform = TRUE,
                             plot = TRUE,
                             focus_color = "skyblue",
                             bottleneck_names = TRUE,
                             bottleneck_cex = 3,
                             gg_extra = list()) {

  method <- match.arg(method)

  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
  } else if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be either an igraph object or a data.frame representing an edge list.")
  }

  if (is.null(vertex_attr(graph, "name"))) {
    vertex_attr(graph, "name") <- as.character(seq_along(1:vcount(graph)))
  }

  # Compute degree and betweenness
  deg <- degree(graph)
  btw <- betweenness(graph, normalized = TRUE)

  # Optional log transform
  deg_val <- if (log_transform) log1p(deg) else deg
  btw_val <- if (log_transform) log1p(btw) else btw

  if (method == "zscore") {
    deg_score <- scale(deg_val)[, 1]
    btw_score <- scale(btw_val)[, 1]
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

  # Compile result
  result <- tibble(
    node = vertex_attr(graph, "name"),
    degree = deg,
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
      list(
        plot = p,
        method = paste("Bottlenecks identified by method:", method,
                       "with Degree metric threshold =", degree_threshold,
                       "and Betweenness metric threshold =", betweenness_threshold),
        result = result,
        graph = graph
      ))

}
