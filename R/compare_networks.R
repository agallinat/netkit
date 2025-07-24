#' Compare Two Networks
#'
#' This function compares two networks using summary metrics, degree distributions,
#' and topological similarity measures. It overlays the complementary cumulative
#' frequency distributions (CCDFs) of degree and returns a combined report.
#'
#' @param graph1 An igraph object or data.frame (edge list).
#' @param graph2 An igraph object or data.frame (edge list).
#' @param remove_singles Logical; remove single nodes before analysis.
#' @param label.size Labels' size in the CCDF plot.
#' @param show_PL Logical; whether to fit and display power law exponents.
#' @param PL_exponents Vector; power-law slopes to show.
#' @param colors Optional vector of colors for the CCDF plot.
#'
#' @return A list with:
#'   - metrics1: Summary metrics for graph1
#'   - metrics2: Summary metrics for graph2
#'   - CCDF_plot: ggplot overlay of ccdfs
#'   - jaccard_similarity: Jaccard index of edge sets
#'   - node_overlap: Fraction of shared nodes
#'   - edge_overlap: Fraction of shared edges
#'   - ks_test: KS test result for degree distributions
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g1 <- sample_pa(100)
#' g2 <- sample_gnp(100, 0.05, directed = F)
#' compare_networks(g1, g2)
#' }
#'
#' @importFrom dplyr bind_rows
#' @importFrom igraph is_igraph degree induced_subgraph graph_from_data_frame V E
#' @importFrom ggplot2 ggplot aes geom_line aes_string scale_color_manual labs coord_cartesian theme_minimal scale_y_log10
#' @importFrom scales trans_breaks trans_format math_format label_math hue_pal
#' @importFrom stats setNames ks.test
#' @importFrom graphics par
#'
#' @export
compare_networks <- function(graph1, graph2,
                             remove_singles = FALSE,
                             show_PL = TRUE,
                             PL_exponents = c(2, 3),
                             colors = c("#e41a1c", "#000831", "#9c52f2", "#b8b8ff"),
                             label.size = 12
) {

  # --- Validate input ---
  if (inherits(graph1, "data.frame")) {
    graph1 <- igraph::graph_from_data_frame(graph1, directed = FALSE)
  } else if (!igraph::is_igraph(graph1)) {
    stop("Input 'graph1' must be an igraph object or edge list data.frame")
  }
  if (inherits(graph2, "data.frame")) {
    graph2 <- igraph::graph_from_data_frame(graph2, directed = FALSE)
  } else if (!igraph::is_igraph(graph2)) {
    stop("Input 'graph2' must be an igraph object or edge list data.frame")
  }

  # --- Handle single nodes ---
  if (remove_singles) {
    degree_g1 <- igraph::degree(graph1)
    graph1 <- induced_subgraph(graph1, vids = which(degree_g1 > 0))
    degree_g2 <- igraph::degree(graph2)
    graph2 <- induced_subgraph(graph2, vids = which(degree_g2 > 0))
    if (0 %in% c(degree_g1, degree_g2)) {
      cat("Single nodes excluded from the analysis.\nSet 'remove_singles' to FALSE to include all nodes.\n")
    }
  }

  # --- Compute basic metrics ---
  metrics1 <- summarize_graph_metrics(graph1)
  metrics2 <- summarize_graph_metrics(graph2)

  # --- Degree distributions ---
  deg1 <- igraph::degree(graph1)
  deg2 <- igraph::degree(graph2)

  # --- KS Test ---
  ks <- suppressWarnings(stats::ks.test(deg1, deg2))

  # --- Similarity ---
  nodes1 <- igraph::vertex_attr(graph1, "name")
  nodes2 <- igraph::vertex_attr(graph2, "name")
  shared_nodes <- intersect(nodes1, nodes2)
  node_overlap <- length(shared_nodes) / length(union(nodes1, nodes2))

  edge_set1 <- apply(igraph::as_edgelist(graph1), 1, function(x) paste(sort(x), collapse = "|"))
  edge_set2 <- apply(igraph::as_edgelist(graph2), 1, function(x) paste(sort(x), collapse = "|"))
  edge_overlap <- length(intersect(edge_set1, edge_set2)) / length(union(edge_set1, edge_set2))
  jaccard_sim <- length(intersect(edge_set1, edge_set2)) / length(union(edge_set1, edge_set2))

  similarity <- data.frame(jaccard_similarity = jaccard_sim,
                           node_overlap = node_overlap,
                           edge_overlap = edge_overlap)

  # Compute ccdfs
  ccdf1 <- compute_ccdf(graph1)
  ccdf1$graph <- "Graph 1"

  ccdf2 <- compute_ccdf(graph2)
  ccdf2$graph <- "Graph 2"

  ccdf_combined <- rbind(ccdf1, ccdf2)

  # Ensure 'colors' has enough entries
  n_lines <- if (show_PL) length(PL_exponents) + 2 else 2

  if (length(colors) < n_lines) {
    colors <- hue_pal(l = 65, c = 100)(n_lines)
  }

  colors_vec <- c("Graph 1" = colors[1], "Graph 2" = colors[2])

  if (show_PL) {
    for (i in seq_along(PL_exponents)) {
      gamma <- PL_exponents[i]
      ccdf_combined[[paste0("PL", gamma)]] <- ccdf_combined$degree^(-gamma)
      colors_vec <- c(colors_vec, setNames(colors[i + 2], paste0("gamma = ", gamma)))
    }
  }

  # Plot overlay
  p_combined <- ggplot(ccdf_combined, aes(x = degree, y = ccdf, color = graph)) +
    geom_line(aes(y=ccdf), linewidth = 0.7)

  if (show_PL) {
    for (gamma in PL_exponents) {
      col_name <- paste0("PL", gamma)
      label <- paste0("gamma = ", gamma)
      p_combined <- p_combined + geom_line(aes_string(y = col_name, color = shQuote(label)),
                                           linetype = "dashed", linewidth = 0.5)
    }
  }

  p_combined <- p_combined +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^floor(x)),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = "log2") +
    labs(x = "Degree, k", y = "Pr(K > k)", color = "Graph") +
    scale_color_manual(values = colors_vec) +
    theme_minimal(base_size = label.size)

  return(list(
    CCDF_plot = p_combined,
    global_topology = rbind(metrics1, metrics2),
    similarity = similarity,
    ks_test = ks
  ))
}
