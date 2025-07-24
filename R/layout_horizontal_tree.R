#' Horizontal Tree Layout for Graph Visualization
#'
#' Rotates a tree layout by -90 degrees to produce a horizontal orientation.
#' This is particularly useful for hierarchical visualizations where a left-to-right structure
#' is preferred over the default top-to-bottom tree layout.
#'
#' @param graph An `igraph` object representing the input graph.
#'
#' @return A numeric matrix with 2 columns representing x and y coordinates of each node in the layout.
#' This matrix can be passed to `plot.igraph()` or other plotting functions.
#'
#' @examples
#' \dontrun{
#' g <- igraph::make_tree(10)
#' coords <- layout_horizontal_tree(g)
#' plot_Net(g, layout = coords)
#' }
#'
#' @export
#'
layout_horizontal_tree <- function(graph) {

  angle = -0.5*pi
  RotMat = matrix(c(cos(angle),
                    sin(angle),
                    -sin(angle),
                    cos(angle)), ncol=2)

  coords <- igraph::layout_as_tree(graph) %*% RotMat

  return(coords)

}
