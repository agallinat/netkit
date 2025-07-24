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
#'   \item{\code{method}}{Description of the method and thresholds used.}
#'   \item{\code{result}}{A \code{tibble} with node name, degree, betweenness, transformed metrics, and hub status.}
#'   \item{\code{graph}}{The original graph with a new vertex attribute \code{is_hub}.}
#' }
#' If \code{plot = TRUE}, a scatter plot of degree vs. betweenness is displayed with
#' hub nodes highlighted.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- sample_pa(100)
#' find_hubs(g, method = "quantile", plot = TRUE)
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

  # Compile result
  result <- tibble(
    node = vertex_attr(graph, "name"),
    degree = deg,
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
    list(
      plot = p,
      method = paste("Hub nodes identified by method:", method,
                     "with Degree metric threshold =", degree_threshold,
                     "and Betweenness metric threshold =", betweenness_threshold),
      result = result,
      graph = graph)
    )

}
