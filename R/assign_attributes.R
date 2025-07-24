#' Assign Vertex and Edge Attributes to an igraph Graph
#'
#' Adds or updates vertex and edge attributes in an \code{igraph} object
#' using user-provided metadata tables. Vertex attributes are matched by the first column of \code{nodes_table},
#' and edge attributes are matched using the first two columns of \code{edge_table}, taking graph direction into account.
#' Only matching nodes and edges are updated. Warnings are issued when there are unmatched entries.
#'
#' @param graph An \code{igraph} object or a data frame containing a symbolic edge list.
#' @param nodes_table Optional. A \code{data.frame} whose first column corresponds to vertex names.
#' @param edge_table Optional. A \code{data.frame} whose first two columns correspond to source and target vertices.
#' @param overwrite Logical. If \code{TRUE}, existing attributes are overwritten. If \code{FALSE}, existing attributes are preserved. Default is \code{TRUE}.
#'
#' @return An \code{igraph} object with added or updated attributes.
#' @export
assign_attributes <- function(graph,
                              nodes_table = NULL,
                              edge_table = NULL,
                              overwrite = TRUE) {

  if (inherits(graph, "data.frame")) {
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
  } else if (!igraph::is_igraph(graph)) {
    stop("Input 'graph' must be either an igraph object or a data.frame representing an edge list.")
  }

  # --- Add vertex attributes ---
  if (!is.null(nodes_table)) {
    if (ncol(nodes_table) < 2) stop("nodes_table must have at least two columns.")

    vertex_ids <- as.character(nodes_table[[1]])
    attr_data <- nodes_table[, -1, drop = FALSE]
    graph_nodes <- igraph::V(graph)$name

    # Match nodes_table to graph nodes by name
    match_idx <- match(graph_nodes, vertex_ids)
    n_matched <- sum(!is.na(match_idx))

    if (n_matched < nrow(nodes_table)) {
      warning(sprintf("Only %d of %d nodes in 'nodes_table' matched the graph vertices.", n_matched, nrow(nodes_table)))
    }
    if (n_matched < length(graph_nodes)) {
      warning(sprintf("Only %d of %d graph vertices matched 'nodes_table' entries.", n_matched, length(graph_nodes)))
    }

    existing_v_attrs <- igraph::vertex_attr_names(graph)

    for (col in colnames(attr_data)) {
      if (col %in% existing_v_attrs && !overwrite) next
      if (col %in% existing_v_attrs && overwrite) warning(sprintf("Overwriting existing vertex attribute '%s'", col))

      values <- rep(NA, igraph::vcount(graph))
      values[!is.na(match_idx)] <- attr_data[[col]][match_idx[!is.na(match_idx)]]
      igraph::vertex_attr(graph, col) <- values
    }
  }

  # --- Add edge attributes ---
  if (!is.null(edge_table)) {
    if (ncol(edge_table) < 3) stop("edge_table must have at least three columns (first two = endpoints, others = attributes).")

    from <- as.character(edge_table[[1]])
    to <- as.character(edge_table[[2]])
    attr_data <- edge_table[, -c(1,2), drop = FALSE]

    # Prepare edge data from graph
    edge_ends <- igraph::ends(graph, igraph::E(graph), names = TRUE)
    graph_df <- data.frame(FROM = edge_ends[,1],
                           TO = edge_ends[,2],
                           index = seq_len(nrow(edge_ends)),
                           stringsAsFactors = FALSE)

    # Prepare attribute table
    table_df <- edge_table
    colnames(table_df)[1:2] <- c("FROM", "TO")

    if (!igraph::is_directed(graph)) {
      # Sort endpoints for undirected graphs
      graph_df[, c("FROM", "TO")] <- t(apply(graph_df[, c("FROM", "TO")], 1, sort))
      table_df[, c("FROM", "TO")] <- t(apply(table_df[, c("FROM", "TO")], 1, sort))
    }

    # Merge based on FROM and TO to get matching attributes
    merged <- merge(graph_df, table_df, by = c("FROM", "TO"), all.x = TRUE, sort = FALSE)

    # Report matching stats
    n_matched <- sum(!is.na(merged$index))
    if (n_matched < nrow(table_df)) {
      warning(sprintf("Only %d of %d edges in 'edge_table' matched the graph edges.", n_matched, nrow(table_df)))
    }
    if (igraph::ecount(graph) != n_matched) {
      warning(sprintf("Only %d rows from 'edge_table' matched graph edges (total edges in graph: %d).", n_matched, igraph::ecount(graph)))
    }

    # Assign edge attributes
    existing_e_attrs <- igraph::edge_attr_names(graph)
    for (col in colnames(attr_data)) {
      if (col %in% existing_e_attrs && !overwrite) next
      if (col %in% existing_e_attrs && overwrite) {
        warning(sprintf("Overwriting existing edge attribute '%s'", col))
      }
      values <- rep(NA, igraph::ecount(graph))
      values[merged$index] <- merged[[col]]
      igraph::edge_attr(graph, col) <- values
    }
  }

  return(graph)
}

