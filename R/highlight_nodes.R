#' Highlight Nodes in a Network Plot
#'
#' Highlights selected nodes in an `igraph` object by changing their label, fill color, and/or outline color.
#' The function modifies node attributes and visualizes the result using `plot_Net()`.
#'
#' @param graph An `igraph` object or a `data.frame` representing a symbolic edge list. If a `data.frame`, it should have at least two columns specifying source and target nodes.
#' @param nodes A character vector of node names to highlight.
#' @param method Character vector specifying how to highlight the nodes. One or more of: `"label"`, `"fill"`, `"outline"`.
#' @param label_color Color for highlighted node labels. Used only if `"label"` is in `method`. Default is `"darkred"`.
#' @param highlight_fill_color Fill color for highlighted nodes. Used only if `"fill"` is in `method`. Default is `"orange"`.
#' @param highlight_frame_color Outline color for highlighted nodes. Used only if `"outline"` is in `method`. Default is `"darkblue"`.
#' @param background_color Fill or frame color for non-highlighted nodes. Default is `"gray"`.
#' @param ... Additional arguments passed to `plot_Net()`.
#'
#' @return Invisibly returns the `igraph` object with updated attributes.
#'
#' @examples
#' \dontrun{
#' g <- igraph::make_ring(10)
#' igraph::V(g)$name <- letters[1:10]
#' highlight_nodes(g, nodes = c("a", "j"), method = c("label", "fill"))
#' }
#'
#' @export
highlight_nodes <- function(graph,
                            nodes,
                            method = c("label", "fill", "outline"),
                            label_color = "darkred",
                            highlight_fill_color = "orange",
                            highlight_frame_color = "darkblue",
                            background_color = "gray", ...) {

  # Input checks
  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
  } else if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be either an igraph object or a data.frame representing an edge list.")
  }

  # Ensure node names exist
  if (is.null(igraph::V(graph)$name)) {
    igraph::V(graph)$name <- as.character(seq_len(igraph::vcount(graph)))
  }

  # Logical vector: which nodes to highlight
  highlight_mask <- igraph::V(graph)$name %in% nodes

  if ("label" %in% method) {
    igraph::V(graph)$label <- ifelse(highlight_mask, igraph::V(graph)$name, NA)
  }

  if ("fill" %in% method) {
    igraph::V(graph)$color <- ifelse(highlight_mask, highlight_fill_color, background_color)
  } else {
    igraph::V(graph)$color <- background_color
  }

  if ("outline" %in% method) {
    igraph::V(graph)$frame.color <- ifelse(highlight_mask, highlight_frame_color, background_color)
  } else {
    igraph::V(graph)$frame.color <- NA
  }

  plot_Net(graph,
           color = "color",
           label = "label" %in% method,
           label.color = label_color,
           frame.color = "frame.color", ...)

}
