#' Assign Vertex and Edge Attributes to an igraph Graph
#'
#' Adds or updates vertex and edge attributes in an \code{igraph} object
#' using user-provided metadata tables. Vertex attributes are matched by the first column of \code{nodes_table},
#' and edge attributes are matched using the first two columns of \code{edge_table}, taking graph direction into account.
#' Only matching nodes and edges are updated. Warnings are issued when there are unmatched entries.
#'
#' @section Parallel edges:
#' An endpoint pair does not identify an edge in a multigraph. When a graph holds
#' several edges between the same two vertices, the k-th such edge takes the k-th
#' `edge_table` row carrying that pair, in the order both appear. A graph with no
#' parallel edges is unaffected, since every pair then occurs once on each side.
#' Rows and edges left over -- a pair repeated more often on one side than the
#' other -- are reported by the usual unmatched warnings, and unmatched edges
#' receive `NA`.
#'
#' In an undirected graph the endpoint pair is read as unordered, so an
#' `edge_table` row may list either endpoint first.
#'
#' @param graph An \code{igraph} object or a data frame containing a symbolic edge list.
#' @param nodes_table Optional. A \code{data.frame} whose first column corresponds to vertex names.
#' @param edge_table Optional. A \code{data.frame} whose first two columns correspond to source and target vertices.
#' @param overwrite Logical. If \code{TRUE}, existing attributes are overwritten. If \code{FALSE}, existing attributes are preserved. Default is \code{TRUE}.
#'
#' @return An \code{igraph} object with added or updated attributes.
#'
#' @examples
#' g <- igraph::make_ring(5)
#' igraph::V(g)$name <- letters[1:5]
#'
#' nodes <- data.frame(node = letters[1:5], group = c("a", "a", "b", "b", "b"))
#' edges <- data.frame(from = c("a", "b"), to = c("b", "c"), weight = c(10, 20))
#'
#' g <- assign_attributes(g, nodes_table = nodes, edge_table = edges)
#' igraph::vertex_attr(g, "group")
#' igraph::edge_attr(g, "weight")
#'
#' @export
assign_attributes <- function(graph,
                              nodes_table = NULL,
                              edge_table = NULL,
                              overwrite = TRUE) {

  graph <- as_netkit_graph(graph)

  # --- Add vertex attributes ---
  if (!is.null(nodes_table)) {
    if (ncol(nodes_table) < 2) stop("nodes_table must have at least two columns.")

    vertex_ids <- as.character(nodes_table[[1]])
    attr_data <- nodes_table[, -1, drop = FALSE]
    graph_nodes <- igraph::V(graph)$name

    # Match nodes_table to graph nodes by name
    match_idx <- match(graph_nodes, vertex_ids)

    # Both counts are derived from match_idx, so they describe one matching
    # rather than two. Counting matched vertices and then reporting it against
    # nrow(nodes_table) overstates the unmatched rows whenever the table repeats
    # an id, because match() resolves every repeat to the first row.
    n_vertices_matched <- sum(!is.na(match_idx))
    n_rows_used <- length(unique(match_idx[!is.na(match_idx)]))

    if (n_rows_used < nrow(nodes_table)) {
      warning(sprintf("Only %d of %d rows in 'nodes_table' matched a graph vertex.",
                      n_rows_used, nrow(nodes_table)))
    }
    if (n_vertices_matched < length(graph_nodes)) {
      warning(sprintf("Only %d of %d graph vertices matched 'nodes_table' entries.",
                      n_vertices_matched, length(graph_nodes)))
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

    if (igraph::is_directed(graph)) {
      g_a <- edge_ends[, 1]; g_b <- edge_ends[, 2]
      t_a <- from;           t_b <- to
    } else {
      # Endpoint order carries no information in an undirected graph, so both
      # sides are keyed on the sorted pair. pmin/pmax rather than apply(, 1, sort):
      # one vectorized pass instead of one call per edge.
      g_a <- pmin(edge_ends[, 1], edge_ends[, 2])
      g_b <- pmax(edge_ends[, 1], edge_ends[, 2])
      t_a <- pmin(from, to)
      t_b <- pmax(from, to)
    }

    # Parallel edges share an endpoint pair, so the pair alone does not identify
    # an edge. The k-th graph edge of a pair takes the k-th edge_table row of
    # that pair: a 1:1 matching that reduces to plain pair matching when there
    # are no parallel edges. merge() cannot be used here -- it forms the
    # cartesian product of the two duplicate groups, which both inflates the
    # reported counts and, because `values[index] <- ` keeps the last of a
    # duplicated index, assigned every parallel edge whichever row the merge
    # happened to put last.
    #
    # The separator cannot occur in a vertex name read from a table, so the
    # pasted key is unambiguous.
    sep <- "\r"
    g_key <- paste(g_a, g_b, sep = sep)
    t_key <- paste(t_a, t_b, sep = sep)
    g_occ <- stats::ave(seq_along(g_key), g_key, FUN = seq_along)
    t_occ <- stats::ave(seq_along(t_key), t_key, FUN = seq_along)

    row_of_edge <- match(paste(g_key, g_occ, sep = sep),
                         paste(t_key, t_occ, sep = sep))
    matched <- !is.na(row_of_edge)

    # Both counts come from row_of_edge, so a row count and an edge count are
    # never compared against different denominators.
    n_edges_matched <- sum(matched)
    n_rows_used <- length(unique(row_of_edge[matched]))

    if (n_rows_used < nrow(edge_table)) {
      warning(sprintf("Only %d of %d rows in 'edge_table' matched a graph edge.",
                      n_rows_used, nrow(edge_table)))
    }
    if (n_edges_matched < igraph::ecount(graph)) {
      warning(sprintf("Only %d of %d graph edges matched an 'edge_table' row.",
                      n_edges_matched, igraph::ecount(graph)))
    }

    # Assign edge attributes
    existing_e_attrs <- igraph::edge_attr_names(graph)
    for (col in colnames(attr_data)) {
      if (col %in% existing_e_attrs && !overwrite) next
      if (col %in% existing_e_attrs && overwrite) {
        warning(sprintf("Overwriting existing edge attribute '%s'", col))
      }
      values <- rep(NA, igraph::ecount(graph))
      values[matched] <- attr_data[[col]][row_of_edge[matched]]
      igraph::edge_attr(graph, col) <- values
    }
  }

  return(graph)
}

