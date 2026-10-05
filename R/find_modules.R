#' Detect and Visualize Network Modules (Communities)
#'
#' Identifies modules (communities) in a network using a variety of community detection
#' algorithms from the \pkg{igraph} package. Optionally filters out small modules,
#' visualizes the detected modules, and returns induced subgraphs for each module.
#'
#' @param graph An \code{igraph} object representing the network to analyze or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param method Character. Community detection method. Options include:
#'   \code{"louvain"}, \code{"walktrap"}, \code{"infomap"}, \code{"edge_betweenness"},
#'   \code{"fluid_communities"}, \code{"fast_greedy"}, \code{"leading_eigen"},
#'   \code{"leiden"}, and \code{"spinglass"}.
#' @param min_size Integer. Minimum number of nodes required to retain a module.
#'   Modules smaller than this size are discarded. Default is 3.
#' @param no.of.communities Integer. Required only when \code{method = "fluid_communities"}.
#'   Specifies the number of communities to find.
#' @param return_subgraphs Logical. If \code{TRUE}, returns a list of induced subgraphs
#'   for each detected module.
#' @param plot Logical. If \code{TRUE}, generates a network plot colored by module.
#' @param label Logical. If \code{TRUE}, displays node labels in the plot.
#' @param ... Additional parameters passed to the \code{plot_Net()} function for customizing the plot.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{result}}{A tibble mapping each node to its module assignment.}
#'   \item{\code{module_table}}{Deprecated alias for \code{result}, kept for
#'     backward compatibility.}
#'   \item{\code{n_modules}}{The number of modules that meet the \code{min_size} threshold.}
#'   \item{\code{subgraphs}}{A named list of subgraphs for each module (only if \code{return_subgraphs = TRUE}).}
#'   \item{\code{method}}{The community detection method used.}
#'   \item{\code{graph}}{The input graph with assigned 'module' and 'color' as vertex attributes, if \code{plot = TRUE}.}
#' }
#' If \code{plot = TRUE}, a network plot is displayed with nodes colored by module.
#'
#' @details
#' This function is a wrapper around several \pkg{igraph} community detection algorithms,
#' including Louvain (\code{cluster_louvain()}), Walktrap, Infomap, Fast Greedy, and others.
#' It simplifies their application and offers optional filtering, visualization via \code{plot_Net()},
#' and module subgraph extraction.
#'
#' @references
#' Csardi G, Nepusz T. The igraph software package for complex network research.
#' InterJournal, Complex Systems. 2006;1695. \url{https://igraph.org}
#'
#' @examples
#' g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
#' res <- find_modules(g, method = "louvain", plot = FALSE)
#' res$n_modules
#' head(res$module_table)
#'
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. Community detection reads a weight as a
#'   *strength*, except `method = "edge_betweenness"`, which needs costs and is
#'   given them. `method = "fluid_communities"` cannot use weights and warns.
#'   See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @importFrom igraph is_igraph is_directed as_undirected cluster_louvain cluster_walktrap cluster_infomap cluster_edge_betweenness cluster_fluid_communities cluster_fast_greedy cluster_leading_eigen cluster_leiden cluster_spinglass membership induced_subgraph vertex_attr vertex_attr<- layout_with_fr vcount
#' @importFrom dplyr filter %>%
#' @importFrom tibble as_tibble
#' @importFrom scales hue_pal
#' @importFrom utils modifyList
#'
#' @export
find_modules <- function(graph,
                         method = "louvain",
                         min_size = 3,
                         no.of.communities = NULL, # only required for 'fluid_communities' method.
                         return_subgraphs = FALSE,
                         plot = TRUE,
                         label = FALSE,
                         weights = NULL,
                         weight_type = c("strength", "distance"), ...) {

  # --- Validate input ---
  graph <- as_netkit_graph(graph, backfill_names = TRUE)
  w <- as_netkit_weights(graph, weights, weight_type)

  if (w$weighted && method == "fluid_communities") {
    warning("Method 'fluid_communities' does not support edge weights; ",
            "they will be ignored.", call. = FALSE)
  }

  # --- Community detection ---
  #
  # Almost every algorithm here reads a weight as a *strength*: a heavier edge
  # binds its endpoints more tightly into a module. `edge_betweenness` is the
  # exception, because it works by removing high-betweenness edges and therefore
  # needs edge *costs*. Passing the strength vector to it would invert the
  # algorithm's notion of which edges are central -- which is precisely the kind
  # of plausible-looking wrong answer the weight contract exists to prevent.
  if (method == "louvain") {
    if (is_directed(graph)) {
      collapsed <- collapse_to_undirected(graph, w)
      graph <- collapsed$graph
      w <- collapsed$w
      message("Input graph converted to undirected for Louvain clustering.\n")
    }
    comm_result <- cluster_louvain(graph, weights = w$strength)
  } else if (method == "walktrap") {
    comm_result <- cluster_walktrap(graph, weights = w$strength)
  } else if (method == "infomap") {
    # infomap names the argument e.weights rather than weights.
    comm_result <- cluster_infomap(graph, e.weights = w$strength)
  } else if (method == "edge_betweenness") {
    # Costs, not strengths: see the note above.
    comm_result <- cluster_edge_betweenness(graph, weights = w$distance)
  } else if (method == "fluid_communities") {
    # argument "no.of.communities" required
    comm_result <- cluster_fluid_communities(graph, no.of.communities)
  } else if (method == "fast_greedy") {
    # only for undirected
    if (is_directed(graph)) {
      collapsed <- collapse_to_undirected(graph, w)
      graph <- collapsed$graph
      w <- collapsed$w
      message("Input graph converted to undirected for Fast Greedy clustering.")
    }
    comm_result <- cluster_fast_greedy(graph, weights = w$strength)
  } else if (method == "leading_eigen") {
    # only for undirected
    if (is_directed(graph)) {
      collapsed <- collapse_to_undirected(graph, w)
      graph <- collapsed$graph
      w <- collapsed$w
      message("Input graph converted to undirected for Leading Eigen clustering.")
    }
    comm_result <- cluster_leading_eigen(graph, weights = w$strength)
  } else if (method == "leiden") {
    # only for undirected
    if (is_directed(graph)) {
      collapsed <- collapse_to_undirected(graph, w)
      graph <- collapsed$graph
      w <- collapsed$w
      message("Input graph converted to undirected for Leiden clustering.")
    }
    comm_result <- cluster_leiden(graph, weights = w$strength)
  } else if (method == "spinglass") {
    # Extract the Largest Connected Component (LCC)
    comps <- components(graph)
    largest_comp_nodes <- which(comps$membership == which.max(comps$csize))
    # Carry the strengths through the subgraph as an attribute rather than by
    # position: induced_subgraph() drops the edges leaving the component.
    graph_lcc <- induced_subgraph(
      if (w$weighted) igraph::set_edge_attr(graph, ".netkit_s", value = w$strength) else graph,
      largest_comp_nodes
    )

    # Apply spinglass on the LCC
    comm_result <- cluster_spinglass(
      graph_lcc,
      weights = if (w$weighted) igraph::edge_attr(graph_lcc, ".netkit_s") else NULL
    )

    if (length(comps$csize) > 1) {
      warning("Method 'spinglass' cannot work with unconnected graph. Performing analysis on the LCC...")
    }

  }

  # --- Assign modules ---
  mem <- membership(comm_result)
  node_names <- names(mem)

  module_df <- data.frame(
    node = node_names,
    module = as.integer(mem),
    stringsAsFactors = FALSE
  )

  # --- Filter small modules ---
  module_sizes <- table(module_df$module)
  valid_modules <- as.integer(names(module_sizes[module_sizes >= min_size]))
  module_df <- module_df %>% dplyr::filter(module %in% valid_modules)

  # --- Optional subgraphs ---
  subgraph_list <- NULL
  if (return_subgraphs) {
    subgraph_list <- lapply(valid_modules, function(mod_id) {
      induced_subgraph(graph, vids = module_df$node[module_df$module == mod_id])
    })
    names(subgraph_list) <- paste0("module_", valid_modules)
  }

  vertex_attr(graph, "module") <- NA
  vertex_attr(graph, "module")[match(module_df$node, vertex_attr(graph, "name"))] <- module_df$module

  # --- Optional plot ---
  if (plot) {

    module_levels <- sort(unique(na.omit(vertex_attr(graph, "module"))))
    module_map <- setNames(seq_along(module_levels), module_levels)
    mapped_module_ids <- module_map[as.character(vertex_attr(graph, "module"))]

    colors <- rep("gray80", vcount(graph))
    colors[!is.na(mapped_module_ids)] <- scales::hue_pal(l = 65, c = 100)(length(module_levels))[mapped_module_ids[!is.na(mapped_module_ids)]]
    vertex_attr(graph, "color") <- colors

    if (is.null(vertex_attr(graph, "label"))) {
      vertex_attr(graph, "label") <- vertex_attr(graph, "name")
    }

    default_args <- list(graph,
                         label = label,
                         color = "color",
                         node.size.factor = 2,
                         node.degree.map = TRUE,
                         edge.color = "#999999",
                         edge.width.factor = 1,
                         edge.bw.map = TRUE,
                         label.color = 'black',
                         label.size = 1,
                         layout = NULL
    )

    # Merge defaults with any user-provided arguments in ...
    user_args <- list(...)
    plot_args <- modifyList(default_args, user_args)

    # Call plot_Net with combined args
    do.call(plot_Net, plot_args)

  }

  return(as_netkit_result(list(
    # `result` is the package-wide name for the node-level table; `module_table`
    # is kept as an alias for backward compatibility and is deprecated.
    result = as_tibble(module_df),
    module_table = as_tibble(module_df),
    n_modules = length(valid_modules),
    subgraphs = if (return_subgraphs) subgraph_list else NULL,
    method = method,
    graph = graph
  ), "find_modules"))

}
