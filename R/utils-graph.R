#' Validate and normalize a graph argument
#'
#' Internal helper implementing netkit's shared graph-input contract. Nearly
#' every public function opens with the same three steps, which this centralises:
#' accept an `igraph` object or a data.frame edge list, convert the latter, error
#' on anything else, and optionally backfill `V(graph)$name` with vertex indices
#' when the graph has no names.
#'
#' @param graph An `igraph` object, or a data.frame whose first two columns are a
#'   symbolic edge list (further columns become edge attributes).
#' @param directed Passed to [igraph::graph_from_data_frame()] when `graph` is a
#'   data.frame. Ignored when `graph` is already an `igraph` object, whose own
#'   directedness is always preserved.
#' @param arg Name of the caller's argument, interpolated into the error message
#'   so that multi-graph functions can report which argument was wrong.
#' @param backfill_names If `TRUE`, vertices with no `name` attribute are named by
#'   their index. Callers that key results by vertex name need this; callers that
#'   do not should leave it `FALSE`, because adding names is observable (for
#'   example `plot_Net()` falls back to `name` when choosing vertex labels).
#'
#' @return An `igraph` object.
#'
#' @keywords internal
#' @noRd
as_netkit_graph <- function(graph,
                            directed = FALSE,
                            arg = "graph",
                            backfill_names = FALSE) {

  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = directed)
  } else if (!igraph::is_igraph(graph)) {
    stop(
      sprintf(
        "Input '%s' must be either an igraph object or a data.frame representing an edge list.",
        arg
      ),
      call. = FALSE
    )
  }

  if (backfill_names && is.null(igraph::vertex_attr(graph, "name"))) {
    igraph::vertex_attr(graph, "name") <- as.character(seq_len(igraph::vcount(graph)))
  }

  graph
}
