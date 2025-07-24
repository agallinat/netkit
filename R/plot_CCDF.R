#' Plot Complementary Cumulative Degree Distribution (CCDF)
#'
#' This function plots the complementary cumulative distribution function (CCDF)
#' of node degrees in a network and optionally overlays power-law reference curves.
#'
#' @param graph An \code{igraph} object representing the network to analyze or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param keep_direction Logical. Only for directed graphs. If \code{TRUE}, CCDF curves are drawn for 'in'-dgree,
#'   'out'-degree, and 'all'-degree distributions. \code{FALSE} to ignore directionality.
#' @param remove_singles Logical. If \code{TRUE}, nodes with degree 0 are removed from the graph
#'   before computing the CCDF. Default is \code{FALSE}.
#' @param show_PL Logical. If \code{TRUE}, overlays theoretical power-law reference lines of the form
#'   \eqn{P(K > k) \sim k^{-\gamma}}. Default is \code{TRUE}.
#' @param PL_exponents Numeric vector. The \eqn{\gamma} exponents for the power-law curves. Default is \code{c(2, 3)}.
#' @param colors Optional character vector. Custom colors for the graph curve and power-law lines.
#'   If \code{NULL}, default colors are used.
#' @param label.size Numeric. Font size for axis labels and theme. Passed to \code{theme_minimal}.
#'
#' @return A \code{ggplot2} object showing the CCDF of node degrees on a log-log scale.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- sample_pa(1000)
#' plot_CCDF(g, remove_singles = TRUE)
#' }
#'
#' @importFrom igraph is_igraph degree induced_subgraph
#' @importFrom ggplot2 ggplot aes geom_line aes_string scale_color_manual labs coord_cartesian theme_minimal scale_y_log10
#' @importFrom scales trans_breaks trans_format math_format label_math hue_pal
#' @importFrom stats setNames
#' @importFrom graphics par
#' @importFrom rlang sym
#'
#' @export
plot_CCDF <- function(graph,
                      keep_direction = TRUE,
                      remove_singles = FALSE,
                      show_PL = TRUE,
                      PL_exponents = c(2, 3),
                      colors = c("#000831","#e41a1c","darkgreen", "#9c52f2", "#b8b8ff"),
                      label.size = 12) {

  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = keep_direction)
  } else if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be either an igraph object or a data.frame representing an edge list.")
  }

  if (remove_singles) {
    deg_all <- igraph::degree(graph)
    graph <- induced_subgraph(graph, vids = which(deg_all > 0))
    if (0 %in% deg_all) {
      cat("Single nodes excluded from the analysis.\nSet 'remove_singles' to FALSE to include all nodes.\n")
    }
  }

  deg <- igraph::degree(graph, mode = "all")

  result <- compute_ccdf(graph,
                         mode = "all",
                         remove_singles = remove_singles)

  if(0 %in% deg) {
    cat(paste0(round(100*(1-result$ccdf[1]), digits = 4), "% of single nodes find in the network.\n",
               "Set 'remove_singles' to TRUE to exclude them for the analysis.\n"))
  }

  if (keep_direction && is_directed(graph)) {

    result_in <- compute_ccdf(graph,
                              mode = "in",
                              remove_singles = remove_singles)

    colnames(result_in) <- c("degree", "ccdf_in")

    result_out <- compute_ccdf(graph,
                               mode = "out",
                               remove_singles = remove_singles)

    colnames(result_out) <- c("degree", "ccdf_out")

    result <- dplyr::full_join(result, result_in, by = "degree")
    result <- dplyr::full_join(result, result_out, by = "degree")

    result <- result %>%
      dplyr::filter(degree > 0)

  }

  # Handle colors
  n = ifelse(is_directed(graph), 2, 0) * as.numeric(keep_direction) + length(PL_exponents) * as.numeric(show_PL) + 1

  if (is.null(colors) || length(colors) < n) {
    colors <- c("black", hue_pal(l = 65, c = 100)(n-1))
  }

  colors_vec <- c("Graph" = colors[1])

  if (keep_direction && is_directed(graph)) {
    colors_vec <- c(colors_vec, "degree_in" = colors[2], "degree_out" = colors[3])
  }

  if (show_PL) {
    for (i in seq_along(PL_exponents)) {
      gamma <- PL_exponents[i]
      PL_i <- result$degree^(-gamma)
      PL_i[which(is.infinite(PL_i))] <- NA
      result[[paste0("PL", gamma)]] <- PL_i
      if (keep_direction && is_directed(graph)) {
        colors_vec <- c(colors_vec, setNames(colors[i + 3], paste0("gamma = ", gamma)))
      } else {
        colors_vec <- c(colors_vec, setNames(colors[i + 1], paste0("gamma = ", gamma)))
      }
    }
  }

  p <- ggplot(result, aes(x = degree)) +
    geom_line(aes(y = ccdf, color = "Graph"), linewidth = 0.7)

  if (keep_direction && is_directed(graph)) {
    p <- p + geom_line(aes(y = ccdf_in, color = "degree_in"), linewidth = 0.7) +
             geom_line(aes(y = ccdf_out, color = "degree_out"), linewidth = 0.7)
  }

  if (show_PL) {
    for (gamma in PL_exponents) {
      col_name <- paste0("PL", gamma)
      label <- paste0("gamma = ", gamma)
      p <- p + geom_line(aes_string(y = col_name, color = shQuote(label)),
                         linetype = "dashed", linewidth = 0.5)
    }
  }

  # Add base plot layers
  p <- p +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^floor(x)),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(
      trans = "log2"
    )+
    labs(x = "Degree, k", y = "Pr(K > k)", color = "") +
    scale_color_manual(values = colors_vec) +
    coord_cartesian(ylim = c(min(result$ccdf), NA),
                    xlim = c(1, NA)) +
    theme_minimal(base_size = label.size)

  return(p)

}

#' Compute the Complementary Cumulative Distribution Function (CCDF) of Node Degrees
#'
#' Computes the CCDF of node degrees for a given igraph object. The CCDF is useful
#' for visualizing degree distributions, particularly on log-log plots, to identify
#' power-law or heavy-tailed behaviors. Internal helper function for `plot_CCDF()`
#' and `compare_networks()`.
#'
#' @param graph An igraph object representing the graph.
#' @param mode Character string indicating which degree type to compute.
#'   Options are \code{"all"} (default), \code{"in"}, or \code{"out"}.
#'   Only relevant for directed graphs.
#' @param remove_singles Logical. If \code{TRUE}, nodes with degree zero will be removed
#'   before computing the degree distribution.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{degree}{Integer node degree values.}
#'   \item{ccdf}{Complementary cumulative distribution values (P(X ≥ x)).}
#' }
#'
#' @examples
#' \dontrun{
#' g <- igraph::sample_pa(1000, power = 2.5, directed = FALSE)
#' ccdf_data <- compute_ccdf(g)
#' plot(ccdf_data$degree, ccdf_data$ccdf, log = "xy", type = "l")
#' }
#'
#' @keywords internal
#'
compute_ccdf <- function(graph,
                         mode = c("all", "in", "out"),
                         remove_singles = FALSE) {

  mode <- match.arg(mode)

  if (!igraph::is_igraph(graph)) stop("Input must be an igraph object.")

  if (remove_singles) {
    deg_all <- igraph::degree(graph, mode = "all")
    graph <- induced_subgraph(graph, vids = which(deg_all > 0))
  }

  deg <- igraph::degree(graph, mode = mode)
  deg_tab <- table(factor(deg, levels = 0:max(deg)))
  deg_vals <- as.integer(names(deg_tab))
  ccdf_vals <- rev(cumsum(rev(as.numeric(deg_tab)))) / sum(deg_tab)

  result <- data.frame(degree = deg_vals, ccdf = ccdf_vals)
  result <- result[result$degree > 0, ]  # Remove degree = 0
  return(result)
}
