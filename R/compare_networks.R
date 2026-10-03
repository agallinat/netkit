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
#' \describe{
#'   \item{\code{plot}}{A \code{ggplot2} object overlaying both CCDF curves.}
#'   \item{\code{global_topology}}{A data frame with one row per input graph, as
#'     produced by [summarize_graph_metrics()].}
#'   \item{\code{similarity}}{A one-row data frame with three distinct set
#'     statistics: \code{jaccard_similarity}, the Jaccard index of the edge sets
#'     (shared edges over their union); \code{node_overlap}, the Jaccard index of
#'     the vertex sets; and \code{edge_overlap}, the overlap coefficient of the
#'     edge sets (shared edges over the \emph{smaller} of the two edge sets).
#'     The last is the informative companion to Jaccard when the two networks
#'     differ greatly in size -- a small network nested inside a large one scores
#'     near 1 on overlap and near 0 on Jaccard. Each is \code{NaN} where its
#'     denominator is empty.}
#'   \item{\code{ks_test}}{The Kolmogorov-Smirnov test comparing the two degree
#'     distributions, as returned by [stats::ks.test()].}
#' }
#'
#' @examples
#' g1 <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
#' g2 <- igraph::sample_gnp(80, 0.05, directed = FALSE)
#' res <- compare_networks(g1, g2)
#' res$global_topology
#' res$ks_test
#'
#' @importFrom dplyr bind_rows
#' @importFrom igraph is_igraph degree induced_subgraph graph_from_data_frame V E
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual labs coord_cartesian theme_minimal scale_y_log10 scale_x_continuous
#' @importFrom rlang sym
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
  graph1 <- as_netkit_graph(graph1, arg = "graph1")
  graph2 <- as_netkit_graph(graph2, arg = "graph2")

  # --- Handle single nodes ---
  if (remove_singles) {
    degree_g1 <- igraph::degree(graph1)
    graph1 <- induced_subgraph(graph1, vids = which(degree_g1 > 0))
    degree_g2 <- igraph::degree(graph2)
    graph2 <- induced_subgraph(graph2, vids = which(degree_g2 > 0))
    if (0 %in% c(degree_g1, degree_g2)) {
      message("Single nodes excluded from the analysis. ",
              "Set 'remove_singles' to FALSE to include all nodes.")
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
  shared_edges <- length(intersect(edge_set1, edge_set2))

  # Jaccard and the overlap coefficient answer different questions, and
  # `edge_overlap` used to be a duplicate of `jaccard_similarity` -- the identical
  # expression under a second name. Dividing by the smaller edge set instead makes
  # it the overlap (Szymkiewicz-Simpson) coefficient, which is the informative
  # companion when the two networks differ greatly in size: a small curated
  # network nested inside a large screen scores near 1 here and near 0 on Jaccard.
  union_edges <- length(union(edge_set1, edge_set2))
  min_edges <- min(length(unique(edge_set1)), length(unique(edge_set2)))

  jaccard_sim <- if (union_edges == 0) NaN else shared_edges / union_edges
  edge_overlap <- if (min_edges == 0) NaN else shared_edges / min_edges

  similarity <- data.frame(jaccard_similarity = jaccard_sim,
                           node_overlap = node_overlap,
                           edge_overlap = edge_overlap)

  # Compute ccdfs. rep() rather than a scalar: compute_ccdf() returns zero rows
  # for an edgeless graph, and recycling a length-1 value into a zero-row frame
  # is an error ("replacement has 1 row, data has 0").
  ccdf1 <- compute_ccdf(graph1)
  ccdf1$graph <- rep("Graph 1", nrow(ccdf1))

  ccdf2 <- compute_ccdf(graph2)
  ccdf2$graph <- rep("Graph 2", nrow(ccdf2))

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
      # `!!` forces col_name and label at aes() construction time; see plot_CCDF().
      p_combined <- p_combined + geom_line(aes(y = !!sym(col_name), color = !!label),
                                           linetype = "dashed", linewidth = 0.5)
    }
  }

  p_combined <- p_combined +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^floor(x)),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      transform = "log2") +
    labs(x = "Degree, k", y = "Pr(K > k)", color = "Graph") +
    scale_color_manual(values = colors_vec) +
    theme_minimal(base_size = label.size)

  return(list(
    plot = p_combined,
    global_topology = rbind(metrics1, metrics2),
    similarity = similarity,
    ks_test = ks
  ))
}
