#' Plot an igraph network with customizable node sizes and edge widths
#'
#' This function plots an \code{igraph} network graph with options to
#' customize node size based on vertex degree or a numeric vertex attribute,
#' set node and edge colors, label nodes, and adjust edge widths
#' (optionally mapped to edge betweenness centrality).
#'
#' @param graph An \code{igraph} object representing the network to plot or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param label Logical, or character vector. If \code{FALSE}, no node labels
#'   are shown. If \code{TRUE}, node labels are set to vertex names or existing
#'   labels. If a character vector of length equal to the number of vertices,
#'   used as custom node labels.
#' @param color Character. Color for node fill, or existing vertex attribute name
#'   to map node fill color. Default is \code{'#006d77'}.
#' @param color_ramp A character vector of colors (e.g., \code{c("blue", "white", "red")}) used to create a continuous color ramp
#'   for mapping numeric node attributes to colors. Ignored if a fixed color is provided or if node color is categorical.
#'   Passed to \code{\link[grDevices]{colorRampPalette}} to interpolate a gradient of colors.
#' @param NA_color Node color for NAs in node color mapping.
#' @param frame.color Character. Color for node frame, or existing vertex attribute name
#'   to map node frame color. Default is \code{NULL}.
#' @param node.size.factor Numeric or character. If numeric, acts as a constant
#'   scaling factor for node size. If character, it should be the name of a
#'   numeric vertex attribute used to size nodes.
#' @param node.degree.map Logical. If \code{TRUE} and \code{node.size.factor} is numeric,
#'   node size is mapped to vertex degree multiplied by \code{node.size.factor}.
#'   Ignored if \code{node.size.factor} is character.
#' @param edge.color Character. Color of edges. Default is \code{'#999999'}.
#' @param edge.width.factor Numeric. Factor to scale edge widths. Default is 1.
#' @param edge.bw.map Logical. If \code{TRUE}, edge widths are mapped to the
#'   log-transformed edge betweenness centrality. Otherwise, edges have uniform width.
#' @param label.color Character. Color of node labels. Default is \code{'#e29578'}.
#' @param label.size Numeric. Character expansion factor for label size. Default is 1.
#' @param layout Optional numeric matrix specifying vertex coordinates for layout.
#' @param ... Additional parameters passed to \code{plot.igraph}.
#'
#' @return Invisibly returns \code{NULL}. The function produces a plot.
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- random.graph.game(100, 0.02)
#' plot_Net(g, label = TRUE, node.size.factor = 2)
#' }
#'
#' @importFrom igraph vertex_attr_names vertex_attr degree edge_betweenness vcount vertex_attr<- V<- is_igraph E<- edge_attr edge_attr<-
#' @importFrom graphics par rect text
#' @importFrom stats na.omit
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
plot_Net <- function(
    graph,
    label = FALSE,
    color = '#006d77',
    color_ramp = c("blue","white","red"),
    NA_color = "gray",
    frame.color = NULL,
    node.size.factor = 1,
    node.degree.map = TRUE,
    edge.color = "#999999",
    edge.width.factor = 1,
    edge.bw.map = TRUE,
    label.color = '#e29578',
    label.size = 1,
    layout = NULL, ...
) {

  # Input checks
  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
  } else if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be either an igraph object or a data.frame representing an edge list.")
  }
  if (!is.null(layout) && !is.matrix(layout)) stop("'layout' must be a matrix.")

  # Normalizing function
  safe_normalize <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) {
      return(rep(0.5, length(x)))  # All values same → return neutral size
    } else {
      return((x - rng[1]) / diff(rng))
    }
  }

  # ---- Determine node size factor ----
  if (is.character(node.size.factor)) {
    if (!node.size.factor %in% vertex_attr_names(graph)) {
      stop(paste0("Vertex attribute '", node.size.factor, "' not found in graph."))
    }
    size_values <- vertex_attr(graph, node.size.factor)
    if (!is.numeric(size_values)) stop("Vertex attribute used for node size must be numeric.")
    size_values <- safe_normalize(size_values)
    vertex_attr(graph, "size") <- (size_values + 0.1) * 5  # Rescale like before
    node.degree.map <- FALSE  # Override degree mapping
  }

  # ---- Node size by degree or constant ----
  if (node.degree.map) {
    deg <- degree(graph)
    vertex_attr(graph, "size") <- safe_normalize(deg)
    vertex_attr(graph, "size") <- (vertex_attr(graph, "size") + 0.1) * 5 * node.size.factor
  } else if (!is.character(node.size.factor)) {
    vertex_attr(graph, "size") <- 5 * node.size.factor
  }

  # ---- Vertex aesthetics ----
  # ---- Determine node color factor ----
  legend <- FALSE

  if (color %in% vertex_attr_names(graph)) {
    if (is.numeric(vertex_attr(graph, color))) {

      legend <- TRUE

      vals <- vertex_attr(graph, color)
      vals_norm <- safe_normalize(vals)

      # Create a color palette function from user-defined colors
      color_palette = colorRampPalette(color_ramp)(100)
      pal_fun <- scales::col_numeric(
        palette = color_palette,
        domain = NULL
      )

      # Map normalized values to colors
      V(graph)$color <- ifelse(is.na(vals_norm), NA_color, pal_fun(vals_norm))

    } else {
      vertex_attr(graph, "color") <- vertex_attr(graph, color)
    }
  } else {
    vertex_attr(graph, "color") <- color
  }

  # ---- Outline ----
  if (!is.null(frame.color)) {
    if (color %in% vertex_attr_names(graph)) {
      vertex_attr(graph, "frame.color") <- vertex_attr(graph, frame.color)
    } else {
      vertex_attr(graph, "frame.color") <- frame.color
    }
  } else {
    vertex_attr(graph, "frame.color") <- V(graph)$color
  }

  # ---- Labels ----
  if (is.logical(label)) {
    if (label) {
      if (is.null(vertex_attr(graph, "label")) && !is.null(vertex_attr(graph, "name"))) {
        vertex_attr(graph, "label") <- as.character(vertex_attr(graph, "name"))
        cat("Vertex attribute 'label' is missing, using 'name' for labels.\n")
      } else {
        if (is.null(vertex_attr(graph, "label")) && is.null(vertex_attr(graph, "name"))) {
          vertex_attr(graph, "label") <- as.character(seq_along(1:vcount(graph)))
          cat("Vertex attributes 'label' and 'name' are missing, using vertex index for labels.\n")
        }
      }
    } else {
      vertex_attr(graph, "label") <- NA
    }
  } else if (is.character(label) && length(label) == vcount(graph)) {
    vertex_attr(graph, "label") <- label
  } else {
    stop("Invalid 'label' argument. Must be FALSE, TRUE, or a character vector of same length as number of nodes.")
  }

  vertex_attr(graph, "label.color") <- label.color

  # ---- Edge aesthetics ----
  edge_attr(graph, "color") <- edge.color

  if (edge.bw.map) {
    bw <- log(edge_betweenness(graph, directed = FALSE) + 1) + 1
    bw <- bw / max(bw)
    edge_attr(graph, "width") <- bw * edge.width.factor
  } else {
    edge_attr(graph, "width") <- edge.width.factor
  }

  edge_attr(graph, "weight") <- edge_attr(graph, "width")  # optional

  # ---- Plotting ----
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  # ---- Legend ----
  if(legend) {

    par(mar = c(0, 0, 0, 1))

    plot(
      graph,
      layout = layout,
      vertex.label.cex = label.size, ...
    )

    # Add color legend (simplified: no title, only top/bottom labels)
    legend_gradient <- function(xl, yb, xr, yt, col = color_palette,
                                values = values) {
      par(xpd = TRUE)
      rect_seq <- seq(yb, yt, length.out = length(col) + 1)
      for (i in seq_along(col)) {
        rect(xl, rect_seq[i], xr, rect_seq[i + 1], col = col[i], border = NA)
      }
      # Show only bottom and top labels
      text(x = xr + 0.01, y = c(yb, yt), labels = values, adj = 0)
    }


    # Draw legend on the side (tweak coordinates as needed)
    legend_gradient(xl = 1.1, yb = 0.2, xr = 1.15, yt = 0.8, col = color_palette,
                    values = range(na.omit(vals)))

  } else {

    par(mar = c(0, 0, 0, 0))

    plot(
      graph,
      layout = layout,
      vertex.label.cex = label.size, ...
    )
  }
}
