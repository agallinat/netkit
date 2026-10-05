#' Extract an Interpretable Subnetwork Around a Set of Nodes
#'
#' Given a set of nodes of interest -- a gene list, a set of disease genes, a
#' group of drug targets -- returns a connected, readable subnetwork around them.
#' Five strategies are offered, from the plain induced subgraph through to an
#' approximate Steiner tree and a diffusion-based expansion.
#'
#' @param graph An `igraph` object representing the network to analyze, or a data
#'   frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param nodes Character vector of vertex names to build the subnetwork around
#'   (the "seeds", or in the Steiner case the "terminals").
#' @param method How to choose the nodes to keep:
#'   \describe{
#'     \item{`"induced"`}{Only `nodes` themselves, with whatever edges run between
#'       them. The baseline every other method should be compared against.}
#'     \item{`"neighbors"`}{`nodes` plus their `order`-step neighborhood.}
#'     \item{`"shortest_paths"`}{The union of all shortest paths between every
#'       pair of `nodes`. The classic "connect my gene list".}
#'     \item{`"steiner"`}{An approximate minimum Steiner tree over `nodes`: the
#'       smallest tree that connects them all. Usually the most readable figure,
#'       because it is a tree rather than a union of overlapping paths.}
#'     \item{`"diffusion"`}{Propagate from `nodes` with [network_diffusion()] and
#'       keep the `top_n` highest-scoring vertices.}
#'   }
#' @param order Integer. Neighborhood radius when `method = "neighbors"`.
#'   Default is `1`.
#' @param max_degree Integer or `NULL`. When `method = "neighbors"`, exclude
#'   neighbors whose degree exceeds this. Without it, first-neighbor expansion
#'   on a hub-dominated network returns most of the network; see Details.
#' @param top_n Integer. Number of vertices to keep when `method = "diffusion"`.
#'   Default is `100`.
#' @param diffusion_method Passed to [network_diffusion()] as `method` when
#'   `method = "diffusion"`. Default is `"rwr"`.
#' @param largest_component Logical. If `TRUE`, return only the largest connected
#'   component of the result. Default is `FALSE`, which keeps a forest when the
#'   seeds span several components.
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. Path lengths and the Steiner tree use edge *costs*;
#'   diffusion uses *strengths*. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#' @param plot Logical. If `TRUE` (default), draws the extracted subnetwork with
#'   [plot_Net()], seeds highlighted.
#' @param ... Additional arguments passed to [plot_Net()].
#'
#' @return A list with:
#' \describe{
#'   \item{`result`}{A tibble with one row per retained vertex: `node`, `reason`
#'     (why it was kept -- `"seed"`, `"neighbor"`, `"on_path"`, `"steiner"` or
#'     `"diffused"`), `is_seed`, and `score` (the diffusion score, or `NA` for
#'     the other methods).}
#'   \item{`graph`}{The extracted subgraph, with `is_seed` and `reason` attached
#'     as vertex attributes. For `method = "steiner"` this is the Steiner tree
#'     itself -- exactly the tree's edges. For every other method it is the
#'     subgraph *induced* on the selected vertices, so it carries all edges that
#'     run between them, which is more context but may contain cycles.}
#'   \item{`method`}{A human-readable description of what was done.}
#' }
#'   Note there is no `plot` element: like [plot_Net()], [find_modules()] and
#'   [highlight_nodes()], this function renders through base graphics via
#'   `plot.igraph` and so has no plot object to hand back.
#'
#' @details
#' `"shortest_paths"` and `"steiner"` both answer "how are these nodes connected
#' to each other", but differently. The union of shortest paths keeps every
#' equally short route, so it grows quickly and can be dense; the Steiner tree
#' keeps one connecting structure of minimum total cost, so it is always a tree
#' with `vcount - 1` edges and is far easier to read. Use the union when you care
#' about redundancy of connection, the tree when you want a figure.
#'
#' The Steiner tree problem is NP-hard. This implementation uses the
#' Kou-Markowsky-Berman heuristic: build the metric closure over the terminals
#' (their pairwise shortest-path distances), take a minimum spanning tree of that
#' closure, expand each closure edge back into the actual path it stands for, and
#' prune non-terminal leaves. The result is guaranteed within a factor of
#' `2 - 2/|terminals|` of the true optimum, which in practice is close.
#'
#' Because the problem is only solved approximately, the tree returned is one of
#' possibly several near-minimal trees, and which one depends on how ties between
#' equal-cost paths are broken. That tie-breaking is not the same in the weighted
#' and unweighted code paths -- \pkg{igraph} uses breadth-first search when no
#' weights are given and Dijkstra when they are -- so passing `weights` whose
#' values happen to be all equal can return a different tree of slightly
#' different total cost than omitting them. Both satisfy the approximation bound;
#' neither is "the" Steiner tree. Treat the specific vertex set as one valid
#' answer rather than a canonical one.
#'
#' `max_degree` exists because first-neighbor expansion is dominated by hubs. In
#' a protein interaction network a handful of promiscuous proteins are adjacent
#' to a large fraction of the graph, so the one-step neighborhood of almost any
#' gene list is most of the network. Capping degree removes them and leaves a
#' subnetwork whose edges carry information.
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @references
#' Kou, L., Markowsky, G., & Berman, L. (1981). A fast algorithm for Steiner
#' trees. *Acta Informatica*, 15(2), 141-145. \doi{10.1007/BF00288961}
#'
#' @seealso [network_diffusion()] for the scores behind `method = "diffusion"`,
#'   [find_modules()] for unsupervised community structure.
#'
#' @examples
#' g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
#' igraph::V(g)$name <- paste0("n", seq_len(igraph::vcount(g)))
#' seeds <- c("n5", "n20", "n40", "n60")
#'
#' # The minimal tree connecting the seeds: always vcount - 1 edges.
#' st <- extract_subnetwork(g, seeds, method = "steiner", plot = FALSE)
#' igraph::vcount(st$graph)
#' st$result
#'
#' # Every shortest route between them, which is a superset.
#' sp <- extract_subnetwork(g, seeds, method = "shortest_paths", plot = FALSE)
#' igraph::vcount(sp$graph) >= igraph::vcount(st$graph)
#'
#' # `top_n` is small here to keep the example fast on an 80-node graph.
#' df <- extract_subnetwork(g, seeds, method = "diffusion", top_n = 15,
#'                          plot = FALSE)
#' head(df$result)
#'
#' @importFrom igraph induced_subgraph distances ego mst vcount ecount V
#' @importFrom igraph vertex_attr vertex_attr<- degree components
#' @importFrom igraph shortest_paths delete_vertices graph_from_data_frame
#' @importFrom igraph edge_attr edge_attr_names delete_edge_attr set_edge_attr ends E
#' @importFrom tibble tibble
#' @importFrom utils modifyList
#'
#' @export
extract_subnetwork <- function(graph,
                               nodes,
                               method = c("induced", "neighbors",
                                          "shortest_paths", "steiner",
                                          "diffusion"),
                               order = 1,
                               max_degree = NULL,
                               top_n = 100,
                               diffusion_method = c("rwr", "laplacian", "heat"),
                               largest_component = FALSE,
                               weights = NULL,
                               weight_type = c("strength", "distance"),
                               plot = TRUE,
                               ...) {

  method <- match.arg(method)
  diffusion_method <- match.arg(diffusion_method)

  graph <- as_netkit_graph(graph, backfill_names = TRUE)
  w <- as_netkit_weights(graph, weights, weight_type)

  all_names <- igraph::V(graph)$name
  requested <- as.character(nodes)
  seeds <- intersect(requested, all_names)

  if (length(seeds) == 0) {
    stop("None of the 'nodes' are vertices of the graph.", call. = FALSE)
  }
  if (length(seeds) < length(unique(requested))) {
    warning(sprintf(
      "%d of %d 'nodes' are not vertices of the graph and were dropped.",
      length(unique(requested)) - length(seeds), length(unique(requested))
    ), call. = FALSE)
  }

  # Carry the costs on the graph so that every subsetting step below keeps them
  # aligned with the edges by construction rather than by index arithmetic.
  graph_c <- if (w$weighted) {
    igraph::set_edge_attr(graph, ".netkit_d", value = w$distance)
  } else {
    graph
  }
  cost <- function(g) if (w$weighted) igraph::edge_attr(g, ".netkit_d") else NULL

  keep <- switch(
    method,
    induced = stats::setNames(rep("seed", length(seeds)), seeds),
    neighbors = subnetwork_neighbors(graph, seeds, order, max_degree),
    shortest_paths = subnetwork_shortest_paths(graph_c, seeds, cost(graph_c)),
    steiner = subnetwork_steiner(graph_c, seeds, cost(graph_c)),
    diffusion = subnetwork_diffusion(graph, seeds, top_n, diffusion_method,
                                     weights, weight_type)
  )

  scores <- attr(keep, "scores")
  prebuilt <- attr(keep, "graph")
  keep_names <- names(keep)

  # Most methods select a vertex set and the induced subgraph is the right
  # answer. `"steiner"` is different: its whole point is that the result is a
  # *tree*, and inducing on the tree's vertex set puts back every edge the
  # spanning-tree step deliberately dropped -- which silently returned a
  # subgraph with cycles under a method named after a tree.
  sub <- if (!is.null(prebuilt)) {
    prebuilt
  } else {
    igraph::induced_subgraph(graph, keep_names)
  }

  # Vertex attributes carried only to keep weights aligned internally are an
  # implementation detail and must not leak into the returned graph.
  for (ea in intersect(c(".netkit_d", ".netkit_s"), igraph::edge_attr_names(sub))) {
    sub <- igraph::delete_edge_attr(sub, ea)
  }

  if (largest_component && igraph::vcount(sub) > 0) {
    comps <- igraph::components(sub)
    sub <- igraph::induced_subgraph(
      sub, which(comps$membership == which.max(comps$csize))
    )
    keep_names <- igraph::V(sub)$name
    keep <- keep[keep_names]
    if (!is.null(scores)) scores <- scores[keep_names]
  }

  # Order the table by graph order so it lines up with the vertex attributes.
  ord <- igraph::V(sub)$name
  reason <- unname(keep[ord])
  is_seed <- ord %in% seeds
  # Seeds are reported as seeds whatever route kept them, so `reason` reads as
  # "why is this node here" rather than "which code path touched it".
  reason[is_seed] <- "seed"

  igraph::vertex_attr(sub, "is_seed") <- is_seed
  igraph::vertex_attr(sub, "reason") <- reason

  result <- tibble::tibble(
    node = ord,
    reason = reason,
    is_seed = is_seed,
    score = if (is.null(scores)) NA_real_ else unname(scores[ord])
  )

  if (plot && igraph::vcount(sub) > 0) {
    defaults <- list(
      graph = sub,
      color = "is_seed",
      label = TRUE
    )
    do.call(plot_Net, utils::modifyList(defaults, list(...)))
  }

  # No `plot` element: this function renders via base plot.igraph through
  # plot_Net(), exactly as find_modules() and highlight_nodes() do, so there is
  # no plot object to return.
  list(
    result = result,
    graph = sub,
    method = paste0(
      "Subnetwork around ", length(seeds), " node(s) by method '", method, "'",
      switch(method,
             neighbors = paste0(" (order = ", order,
                                if (is.null(max_degree)) "" else
                                  paste0(", max_degree = ", max_degree), ")"),
             diffusion = paste0(" (", diffusion_method, ", top_n = ", top_n, ")"),
             ""),
      "; ", describe_weights(w),
      if (largest_component) "; largest component only" else "",
      "; ", igraph::vcount(sub), " nodes and ", igraph::ecount(sub), " edges"
    )
  )
}

#' Neighborhood expansion around a seed set
#'
#' Internal helper for `extract_subnetwork()`.
#'
#' @param graph An `igraph` object.
#' @param seeds Character vector of seed names.
#' @param order Neighborhood radius.
#' @param max_degree Degree cap for non-seed neighbors, or `NULL`.
#'
#' @return A named character vector of retention reasons, keyed by vertex name.
#'
#' @keywords internal
#' @noRd
subnetwork_neighbors <- function(graph, seeds, order, max_degree) {

  if (!is.numeric(order) || length(order) != 1 || order < 0) {
    stop("'order' must be a single non-negative number.", call. = FALSE)
  }

  nb <- igraph::ego(graph, order = order, nodes = seeds, mode = "all")
  nb_names <- unique(unlist(lapply(nb, function(v) v$name)))

  extra <- setdiff(nb_names, seeds)

  if (!is.null(max_degree)) {
    # Applied only to the added neighbors. A seed is in the set because the
    # caller asked for it, so dropping one for being a hub would silently answer
    # a different question.
    deg <- igraph::degree(graph, v = extra)
    extra <- extra[deg <= max_degree]
  }

  c(stats::setNames(rep("seed", length(seeds)), seeds),
    stats::setNames(rep("neighbor", length(extra)), extra))
}

#' Union of all shortest paths between every pair of seeds
#'
#' Internal helper for `extract_subnetwork()`.
#'
#' @param graph An `igraph` object carrying `.netkit_d` when weighted.
#' @param seeds Character vector of seed names.
#' @param costs Edge costs, or `NULL`.
#'
#' @return A named character vector of retention reasons.
#'
#' @keywords internal
#' @noRd
subnetwork_shortest_paths <- function(graph, seeds, costs) {

  if (length(seeds) < 2) {
    return(stats::setNames(rep("seed", length(seeds)), seeds))
  }

  on_path <- character(0)
  unreachable <- 0L

  for (s in seeds) {
    targets <- setdiff(seeds, s)
    sp <- suppressWarnings(
      igraph::shortest_paths(graph, from = s, to = targets, weights = costs,
                             output = "vpath")
    )
    for (vp in sp$vpath) {
      if (length(vp) == 0) {
        unreachable <- unreachable + 1L
      } else {
        on_path <- c(on_path, vp$name)
      }
    }
  }

  if (unreachable > 0) {
    warning(sprintf(
      paste0("%d seed pair(s) are in different components, so no path connects ",
             "them. The result is a forest rather than a connected subnetwork."),
      unreachable %/% 2L + unreachable %% 2L
    ), call. = FALSE)
  }

  on_path <- unique(on_path)
  extra <- setdiff(on_path, seeds)

  c(stats::setNames(rep("seed", length(seeds)), seeds),
    stats::setNames(rep("on_path", length(extra)), extra))
}

#' Approximate minimum Steiner tree over a terminal set
#'
#' Internal helper for `extract_subnetwork()`, implementing the
#' Kou-Markowsky-Berman heuristic: metric closure over the terminals, minimum
#' spanning tree of the closure, expansion of each closure edge back to its
#' underlying path, then pruning of non-terminal leaves. Guaranteed within
#' `2 - 2/|terminals|` of optimal.
#'
#' The problem is NP-hard, so an exact solution is not available; the point of
#' this method over `"shortest_paths"` is that the result is a *tree*, which is
#' what makes it readable as a figure.
#'
#' @param graph An `igraph` object carrying `.netkit_d` when weighted.
#' @param terminals Character vector of terminal names.
#' @param costs Edge costs, or `NULL`.
#'
#' @return A named character vector of retention reasons.
#'
#' @keywords internal
#' @noRd
subnetwork_steiner <- function(graph, terminals, costs) {

  if (length(terminals) < 2) {
    out <- stats::setNames(rep("seed", length(terminals)), terminals)
    attr(out, "graph") <- igraph::induced_subgraph(graph, terminals)
    return(out)
  }

  # 1. Metric closure: pairwise shortest-path distances between terminals.
  D <- igraph::distances(graph, v = terminals, to = terminals, weights = costs)

  # Terminals in different components have infinite distance. Solve each
  # component separately and return a forest, which is the honest answer; the
  # alternative is to drop terminals the caller asked for.
  reachable <- is.finite(D)
  if (!all(reachable)) {
    warning(sprintf(
      paste0("Terminals span %d components, so no single tree connects them. ",
             "Returning a Steiner forest -- one tree per component."),
      igraph::components(igraph::induced_subgraph(graph, terminals))$no
    ), call. = FALSE)
  }

  # 2. MST of the closure, computed per component of the closure graph.
  #    Infinite distances are dropped rather than passed to mst(), which cannot
  #    represent them.
  n_t <- length(terminals)
  edges <- list()
  for (i in seq_len(n_t - 1)) {
    for (j in seq(i + 1, n_t)) {
      if (is.finite(D[i, j])) {
        edges[[length(edges) + 1]] <- data.frame(
          from = terminals[i], to = terminals[j], weight = D[i, j],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(edges) == 0) {
    # No two terminals are connected at all: the forest is the terminals alone.
    out <- stats::setNames(rep("seed", n_t), terminals)
    attr(out, "graph") <- igraph::induced_subgraph(graph, terminals)
    return(out)
  }

  closure <- igraph::graph_from_data_frame(do.call(rbind, edges),
                                           directed = FALSE,
                                           vertices = data.frame(name = terminals))
  closure_mst <- igraph::mst(closure, weights = igraph::E(closure)$weight)

  # 3. Expand each MST edge back into the path it stands for.
  on_tree <- terminals
  mst_ends <- igraph::ends(closure_mst, igraph::E(closure_mst), names = TRUE)
  for (k in seq_len(nrow(mst_ends))) {
    sp <- igraph::shortest_paths(graph, from = mst_ends[k, 1],
                                 to = mst_ends[k, 2], weights = costs,
                                 output = "vpath")
    on_tree <- c(on_tree, sp$vpath[[1]]$name)
  }
  on_tree <- unique(on_tree)

  # 4. Take an MST of the induced subgraph, then prune non-terminal leaves
  #    repeatedly. Step 3 can pick up vertices that two expanded paths share, and
  #    the union of paths need not be a tree; the MST makes it one and the
  #    pruning removes the branches that lead nowhere.
  sub <- igraph::induced_subgraph(graph, on_tree)
  sub <- igraph::mst(sub, weights = if (is.null(costs)) NULL else
    igraph::edge_attr(sub, ".netkit_d"))

  repeat {
    deg <- igraph::degree(sub)
    leaves <- igraph::V(sub)$name[deg <= 1]
    drop <- setdiff(leaves, terminals)
    if (length(drop) == 0) break
    sub <- igraph::delete_vertices(sub, drop)
  }

  kept <- igraph::V(sub)$name
  extra <- setdiff(kept, terminals)

  out <- c(stats::setNames(rep("seed", length(terminals)), terminals),
           stats::setNames(rep("steiner", length(extra)), extra))

  # Hand back the tree itself, not merely its vertex names. The caller would
  # otherwise take the induced subgraph and restore the cycle-closing edges that
  # mst() removed, so the "tree" would not be one.
  attr(out, "graph") <- sub
  out
}

#' Diffusion-based subnetwork expansion
#'
#' Internal helper for `extract_subnetwork()`. Reuses the package's existing
#' diffusion machinery to rank vertices by their proximity to the seed set, then
#' keeps the top `top_n`.
#'
#' @param graph An `igraph` object.
#' @param seeds Character vector of seed names.
#' @param top_n Number of vertices to retain.
#' @param diffusion_method Passed to `network_diffusion()`.
#' @param weights,weight_type Passed to `network_diffusion()`.
#'
#' @return A named character vector of retention reasons, with a `scores`
#'   attribute.
#'
#' @keywords internal
#' @noRd
subnetwork_diffusion <- function(graph, seeds, top_n, diffusion_method,
                                 weights, weight_type) {

  if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) {
    stop("'top_n' must be a single number of at least 1.", call. = FALSE)
  }

  scored <- network_diffusion(graph, seed_nodes = seeds,
                              method = diffusion_method,
                              weights = weights, weight_type = weight_type)

  # network_diffusion() returns rows sorted by descending score. The seeds are
  # kept whether or not they rank in the top_n, because the caller asked for
  # them; top_n therefore bounds the *added* vertices plus the seeds.
  n_keep <- min(top_n, nrow(scored))
  top <- scored$node[seq_len(n_keep)]
  kept <- union(seeds, top)
  extra <- setdiff(kept, seeds)

  out <- c(stats::setNames(rep("seed", length(seeds)), seeds),
           stats::setNames(rep("diffused", length(extra)), extra))
  attr(out, "scores") <- stats::setNames(scored$score, scored$node)[names(out)]
  out
}
