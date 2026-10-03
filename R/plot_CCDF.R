#' Plot Complementary Cumulative Degree Distribution (CCDF)
#'
#' This function plots the complementary cumulative distribution function (CCDF)
#' of node degrees in a network and optionally overlays power-law reference curves.
#'
#' @param graph An \code{igraph} object representing the network to analyze or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param keep_direction Logical. Only for directed graphs. If \code{TRUE}, CCDF curves are drawn for 'in'-degree,
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
#' g <- igraph::sample_pa(200, power = 1.5, directed = FALSE)
#' plot_CCDF(g)
#'
#' # Compare against reference power-law slopes.
#' plot_CCDF(g, PL_exponents = c(2, 2.5, 3))
#'
#' @param weights Optional edge weights: `NULL` (default) to plot the degree
#'   distribution, or the name of an edge attribute / a numeric vector of length
#'   `igraph::ecount(graph)` to plot the vertex *strength* distribution instead.
#'   See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @importFrom igraph is_igraph degree induced_subgraph strength
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual labs coord_cartesian theme_minimal scale_y_log10 scale_x_continuous
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
                      label.size = 12,
                      weights = NULL,
                      weight_type = c("strength", "distance")) {

  # An edge-list data.frame is read as directed only when the caller asked to keep
  # direction; an igraph input keeps its own directedness either way.
  graph <- as_netkit_graph(graph, directed = keep_direction)

  if (remove_singles) {
    deg_all <- igraph::degree(graph)
    graph <- induced_subgraph(graph, vids = which(deg_all > 0))
    if (0 %in% deg_all) {
      cat("Single nodes excluded from the analysis.\nSet 'remove_singles' to FALSE to include all nodes.\n")
    }
  }

  deg <- igraph::degree(graph, mode = "all")

  # Attach the resolved strengths to the graph so every compute_ccdf() call below
  # -- including the in/out ones, which subset -- reads the same aligned values.
  w <- as_netkit_weights(graph, weights, weight_type)
  weight_attr <- NULL
  if (w$weighted) {
    weight_attr <- ".netkit_s"
    graph <- igraph::set_edge_attr(graph, weight_attr, value = w$strength)
  }

  result <- compute_ccdf(graph,
                         mode = "all",
                         remove_singles = remove_singles,
                         weight_attr = weight_attr)

  if(0 %in% deg) {
    cat(paste0(round(100*(1-result$ccdf[1]), digits = 4), "% of single nodes find in the network.\n",
               "Set 'remove_singles' to TRUE to exclude them for the analysis.\n"))
  }

  if (keep_direction && is_directed(graph)) {

    result_in <- compute_ccdf(graph,
                              mode = "in",
                              remove_singles = remove_singles,
                              weight_attr = weight_attr)

    colnames(result_in) <- c("degree", "ccdf_in")

    result_out <- compute_ccdf(graph,
                               mode = "out",
                               remove_singles = remove_singles,
                               weight_attr = weight_attr)

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
      # `!!` forces col_name and label at aes() construction time. Without it each
      # layer would capture the loop variables by reference and every reference
      # line would end up using the final exponent.
      p <- p + geom_line(aes(y = !!sym(col_name), color = !!label),
                         linetype = "dashed", linewidth = 0.5)
    }
  }

  # Add base plot layers
  p <- p +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^floor(x)),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(
      transform = "log2"
    )+
    labs(x = if (w$weighted) "Strength, s" else "Degree, k",
         y = if (w$weighted) "Pr(S > s)" else "Pr(K > k)", color = "") +
    scale_color_manual(values = colors_vec) +
    # An edgeless graph yields no positive-degree classes, so `result` is empty and
    # min() would be Inf (with a warning). Leave the limit to ggplot2 in that case.
    coord_cartesian(ylim = c(if (nrow(result) > 0) min(result$ccdf) else NA, NA),
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
#' @keywords internal
#'
compute_ccdf <- function(graph,
                         mode = c("all", "in", "out"),
                         remove_singles = FALSE,
                         weight_attr = NULL) {

  mode <- match.arg(mode)

  if (!igraph::is_igraph(graph)) stop("Input must be an igraph object.")

  if (remove_singles) {
    deg_all <- igraph::degree(graph, mode = "all")
    # Weights ride along as an edge attribute precisely so that this subsetting
    # cannot desynchronise them from the edges.
    graph <- induced_subgraph(graph, vids = which(deg_all > 0))
  }

  # A strength distribution is continuous, so it cannot use the integer-degree
  # tabulation below: factor(levels = 0:max) would need one level per distinct
  # value and max() is not an integer. Evaluate the CCDF on the observed values
  # instead, which is the same definition without the binning assumption.
  if (!is.null(weight_attr)) {
    s_vals <- igraph::strength(graph, mode = mode,
                               weights = igraph::edge_attr(graph, weight_attr))
    n_v <- length(s_vals)
    if (n_v == 0) {
      return(data.frame(degree = numeric(0), ccdf = numeric(0)))
    }
    vals <- sort(unique(s_vals[s_vals > 0]))
    if (length(vals) == 0) {
      return(data.frame(degree = numeric(0), ccdf = numeric(0)))
    }
    # Divided by the same n_v the values came from, matching the unweighted path,
    # which divides by sum(deg_tab) rather than by a separately derived total.
    ccdf_vals <- vapply(vals, function(v) sum(s_vals >= v) / n_v, numeric(1))
    return(data.frame(degree = vals, ccdf = ccdf_vals))
  }

  deg <- igraph::degree(graph, mode = mode)

  # A graph with no vertices, or no edges, has no positive-degree classes and so
  # no degree distribution to describe. Returning the empty table rather than
  # erroring matches summarize_graph_metrics(), which reports NaN on degenerate
  # input instead of refusing it. max() of an empty vector would be -Inf and make
  # the factor levels invalid, so the length guard comes first.
  if (length(deg) == 0 || max(deg) == 0) {
    return(data.frame(degree = integer(0), ccdf = numeric(0)))
  }

  deg_tab <- table(factor(deg, levels = 0:max(deg)))
  deg_vals <- as.integer(names(deg_tab))
  ccdf_vals <- rev(cumsum(rev(as.numeric(deg_tab)))) / sum(deg_tab)

  result <- data.frame(degree = deg_vals, ccdf = ccdf_vals)
  result <- result[result$degree > 0, ]  # Remove degree = 0
  return(result)
}
