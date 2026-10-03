#' Validate and normalize an edge-weight argument
#'
#' Internal helper implementing netkit's shared edge-weight contract, the
#' companion to `as_netkit_graph()`. Functions that can use edge weights call it
#' immediately after validating the graph.
#'
#' The problem it solves is that \pkg{igraph} gives the `weight` edge attribute
#' two incompatible meanings depending on the function, and picks it up
#' implicitly. For path-based measures ([igraph::betweenness()],
#' [igraph::distances()], [igraph::diameter()], [igraph::mean_distance()]) a
#' weight is a *cost*: larger means further apart. For community detection
#' ([igraph::cluster_louvain()] and friends) it is a *strength*: larger means
#' more tightly bound. netkit's own matrix code, meanwhile, built the adjacency
#' matrix with `as_adjacency_matrix()` and no `attr`, so it ignored weights
#' altogether.
#'
#' The practical consequence was that attaching a confidence score -- the most
#' common edge attribute in the domain netkit targets -- silently inverted every
#' path-based metric (a high-confidence edge was treated as a long distance),
#' was read correctly by community detection, and was dropped entirely by
#' diffusion. Three semantics, none documented.
#'
#' This function makes the choice explicit and derives both vectors once, so each
#' metric can be handed the one it actually needs.
#'
#' @param graph An `igraph` object, already normalized by `as_netkit_graph()`.
#' @param weights `NULL` to ignore edge weights (the default throughout netkit);
#'   the name of an edge attribute; or a numeric vector of length
#'   `igraph::ecount(graph)`.
#' @param type Whether the supplied values are a `"strength"` (larger = more
#'   tightly connected: confidence scores, correlations, read counts) or a
#'   `"distance"` (larger = further apart: costs, dissimilarities). netkit
#'   defaults to `"strength"` because that is what the overwhelming majority of
#'   biological edge annotations are.
#' @param arg Name of the caller's argument, interpolated into error messages.
#' @param warn_unused If `TRUE` and `weights` is `NULL` while the graph *does*
#'   carry a `weight` edge attribute, warn that it is being ignored. This makes
#'   the previously silent behaviour loud at the one moment it matters.
#'
#' @return A list with:
#'   \describe{
#'     \item{`weighted`}{`TRUE` if weights are in use. Call sites branch on this.}
#'     \item{`strength`}{Numeric vector of strengths, or `NULL` when unweighted.
#'       Pass to community detection, to the diffusion adjacency matrix, and
#'       anywhere a larger value should mean a stronger connection.}
#'     \item{`distance`}{Numeric vector of costs, or `NULL` when unweighted. Pass
#'       to [igraph::betweenness()], [igraph::distances()] and the other
#'       path-based measures.}
#'     \item{`type`}{The `type` that was supplied.}
#'     \item{`name`}{The attribute name, or `"<numeric vector>"`, for `method`
#'       strings.}
#'   }
#'
#' @keywords internal
#' @noRd
as_netkit_weights <- function(graph,
                              weights = NULL,
                              type = c("strength", "distance"),
                              arg = "weights",
                              warn_unused = TRUE) {

  type <- match.arg(type)

  if (is.null(weights)) {
    if (warn_unused && "weight" %in% igraph::edge_attr_names(graph)) {
      warning(sprintf(
        paste0("Graph has an edge attribute 'weight' which is being ignored. ",
               "Pass %s = \"weight\" to use it, and set weight_type to declare ",
               "whether it is a 'strength' (larger = more tightly connected) or ",
               "a 'distance' (larger = further apart)."),
        arg
      ), call. = FALSE)
    }
    return(list(weighted = FALSE, strength = NULL, distance = NULL,
                type = type, name = NA_character_))
  }

  # --- Resolve to a numeric vector ---
  if (is.character(weights) && length(weights) == 1) {
    if (!weights %in% igraph::edge_attr_names(graph)) {
      stop(sprintf(
        "Edge attribute '%s' not found. Available edge attributes: %s.",
        weights,
        if (length(igraph::edge_attr_names(graph))) {
          paste(igraph::edge_attr_names(graph), collapse = ", ")
        } else {
          "none"
        }
      ), call. = FALSE)
    }
    name <- weights
    w <- igraph::edge_attr(graph, weights)
  } else if (is.numeric(weights)) {
    name <- "<numeric vector>"
    w <- weights
  } else {
    stop(sprintf(
      paste0("Input '%s' must be NULL, the name of an edge attribute, or a ",
             "numeric vector of length ecount(graph)."),
      arg
    ), call. = FALSE)
  }

  # --- Validate ---
  if (length(w) != igraph::ecount(graph)) {
    stop(sprintf(
      "'%s' has length %d but the graph has %d edges.",
      arg, length(w), igraph::ecount(graph)
    ), call. = FALSE)
  }
  if (!is.numeric(w)) {
    stop(sprintf("'%s' must be numeric.", arg), call. = FALSE)
  }
  if (anyNA(w)) {
    stop(sprintf("'%s' contains missing values.", arg), call. = FALSE)
  }
  if (any(!is.finite(w))) {
    stop(sprintf("'%s' contains non-finite values.", arg), call. = FALSE)
  }
  if (any(w < 0)) {
    stop(sprintf(
      paste0("'%s' contains negative values. Weights must be non-negative; a ",
             "signed attribute (such as activation/inhibition) is not an edge ",
             "weight and cannot be used as one."),
      arg
    ), call. = FALSE)
  }
  # An individual zero is meaningful -- a strength of zero is "no connection",
  # giving an infinite distance, which igraph handles. All zeroes is not: there
  # would be no connection anywhere and every derived quantity is degenerate.
  if (length(w) > 0 && all(w == 0)) {
    stop(sprintf("'%s' is zero for every edge, leaving no connections.", arg),
         call. = FALSE)
  }

  # --- Derive the companion vector ---
  if (type == "strength") {
    strength <- w
    # Reciprocal: the canonical strength-to-cost map. A zero strength becomes an
    # infinite cost, which is exactly right -- igraph treats it as unreachable.
    distance <- ifelse(w == 0, Inf, 1 / w)
  } else {
    distance <- w
    # Reflect rather than invert, so a zero distance (identical endpoints) maps
    # to the largest strength instead of to Inf. Shifting by min() keeps the
    # smallest distance at a positive strength rather than at zero, which would
    # otherwise sever the shortest edge in the graph.
    span <- max(w)
    strength <- span - w + (if (span > 0) min(w[w > 0], span) else 1)
  }

  list(weighted = TRUE, strength = as.numeric(strength),
       distance = as.numeric(distance), type = type, name = name)
}

#' Describe a resolved weight specification for a `method` string
#'
#' Internal. Every netkit function reports what it actually did in its `method`
#' element; this renders the weight half of that sentence consistently.
#'
#' @param w The list returned by `as_netkit_weights()`.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
describe_weights <- function(w) {
  if (!w$weighted) {
    return("unweighted")
  }
  sprintf("weighted by '%s' (interpreted as %s)", w$name, w$type)
}

#' Build a diffusion seed vector
#'
#' Internal. Turns a seed specification into the \eqn{f_0} vector that the
#' diffusion kernels are applied to.
#'
#' The binary form -- 1 for a seed, 0 otherwise -- is the default and matches the
#' original behaviour. `seed_weights` instead starts the diffusion from a
#' continuous signal, which is what the network-propagation literature assumes
#' and what diffusing from a differential-expression result requires: the seeds
#' are not equally important, and their magnitudes carry the evidence.
#'
#' @param all_nodes Character vector of every vertex name, in graph order.
#' @param seed_nodes Character vector of seed names.
#' @param seed_weights `NULL` for a binary seed vector, or a numeric vector of
#'   initial values. May be named (matched to `seed_nodes` by name) or unnamed
#'   (taken as parallel to `seed_nodes`).
#'
#' @return A named numeric vector of length `length(all_nodes)`.
#'
#' @keywords internal
#' @noRd
make_seed_vector <- function(all_nodes, seed_nodes, seed_weights = NULL) {

  if (length(seed_nodes) == 0) {
    stop("'seed_nodes' is empty: diffusion needs at least one seed.",
         call. = FALSE)
  }

  seed_nodes <- as.character(seed_nodes)
  unknown <- setdiff(seed_nodes, all_nodes)
  if (length(unknown) == length(seed_nodes)) {
    stop(sprintf(
      "None of the 'seed_nodes' are vertices of the graph. First few: %s.",
      paste(utils::head(unknown, 5), collapse = ", ")
    ), call. = FALSE)
  }
  if (length(unknown) > 0) {
    warning(sprintf(
      "%d of %d 'seed_nodes' are not vertices of the graph and were dropped: %s%s",
      length(unknown), length(seed_nodes),
      paste(utils::head(unknown, 5), collapse = ", "),
      if (length(unknown) > 5) ", ..." else ""
    ), call. = FALSE)
  }

  f0 <- stats::setNames(numeric(length(all_nodes)), all_nodes)

  if (is.null(seed_weights)) {
    f0[intersect(seed_nodes, all_nodes)] <- 1
    return(f0)
  }

  if (!is.numeric(seed_weights)) {
    stop("'seed_weights' must be numeric.", call. = FALSE)
  }
  if (anyNA(seed_weights) || any(!is.finite(seed_weights))) {
    stop("'seed_weights' contains missing or non-finite values.", call. = FALSE)
  }

  if (!is.null(names(seed_weights))) {
    # Named: match by name, which is the safer form and the one to prefer when
    # the seed list came from a table.
    missing_w <- setdiff(seed_nodes, names(seed_weights))
    if (length(missing_w) > 0) {
      stop(sprintf(
        "'seed_weights' is named but has no entry for %d of the 'seed_nodes': %s%s",
        length(missing_w), paste(utils::head(missing_w, 5), collapse = ", "),
        if (length(missing_w) > 5) ", ..." else ""
      ), call. = FALSE)
    }
    keep <- intersect(seed_nodes, all_nodes)
    f0[keep] <- seed_weights[keep]
  } else {
    if (length(seed_weights) != length(seed_nodes)) {
      stop(sprintf(
        paste0("'seed_weights' has length %d but there are %d 'seed_nodes'. ",
               "Supply one value per seed, or name the vector."),
        length(seed_weights), length(seed_nodes)
      ), call. = FALSE)
    }
    keep <- seed_nodes %in% all_nodes
    f0[seed_nodes[keep]] <- seed_weights[keep]
  }

  f0
}

#' Align seed weights to the retained seeds for permutation testing
#'
#' Internal. `network_diffusion_with_pvalues()` reuses the real seed set's
#' magnitudes on each permuted seed set, so the null tests where the seeds sit
#' rather than what they are worth. A `seed_weights` vector named by the real
#' seed names cannot be matched against permuted names, so it is reduced here to
#' an unnamed vector parallel to the retained seeds.
#'
#' @param seed_weights `NULL`, or a numeric vector as accepted by
#'   `make_seed_vector()`.
#' @param requested_seeds The seed names as the caller supplied them.
#' @param kept_seeds The subset of `requested_seeds` present in the graph, in the
#'   order they will be used.
#'
#' @return `NULL`, or an unnamed numeric vector of length `length(kept_seeds)`.
#'
#' @keywords internal
#' @noRd
resolve_perm_seed_weights <- function(seed_weights, requested_seeds, kept_seeds) {

  if (is.null(seed_weights)) {
    return(NULL)
  }
  if (!is.numeric(seed_weights)) {
    stop("'seed_weights' must be numeric.", call. = FALSE)
  }

  if (!is.null(names(seed_weights))) {
    missing_w <- setdiff(kept_seeds, names(seed_weights))
    if (length(missing_w) > 0) {
      stop(sprintf(
        "'seed_weights' is named but has no entry for %d of the seed nodes: %s%s",
        length(missing_w), paste(utils::head(missing_w, 5), collapse = ", "),
        if (length(missing_w) > 5) ", ..." else ""
      ), call. = FALSE)
    }
    return(unname(seed_weights[kept_seeds]))
  }

  if (length(seed_weights) != length(requested_seeds)) {
    stop(sprintf(
      paste0("'seed_weights' has length %d but there are %d seed nodes. ",
             "Supply one value per seed, or name the vector."),
      length(seed_weights), length(requested_seeds)
    ), call. = FALSE)
  }

  # Drop the entries belonging to seeds that are not in the graph, keeping the
  # correspondence with the retained ones.
  unname(seed_weights[requested_seeds %in% kept_seeds])
}

#' Collapse a directed graph while keeping weights aligned
#'
#' Internal. Several netkit functions coerce a directed graph to undirected with
#' `mode = "collapse"`, which merges reciprocal edge pairs. That renumbers the
#' edges, so an edge-aligned weight vector no longer corresponds to them.
#'
#' Carrying the weights on the graph as edge attributes and letting igraph
#' combine them is the only safe way to do this: a positional index silently
#' misaligns, and the result still looks like a plausible set of weights.
#' Strengths are summed (two reciprocal interactions are stronger than one) and
#' costs are minimised (the cheaper of two routes is the one a path would take).
#'
#' @param graph An `igraph` object.
#' @param w The list returned by `as_netkit_weights()`.
#'
#' @return A list of `graph` (undirected) and `w` (with realigned vectors).
#'
#' @keywords internal
#' @noRd
collapse_to_undirected <- function(graph, w) {

  if (!igraph::is_directed(graph)) {
    return(list(graph = graph, w = w))
  }

  if (!w$weighted) {
    return(list(graph = igraph::as_undirected(graph, mode = "collapse"), w = w))
  }

  graph <- igraph::set_edge_attr(graph, ".netkit_s", value = w$strength)
  graph <- igraph::set_edge_attr(graph, ".netkit_d", value = w$distance)

  # The default edge.attr.comb would drop both attributes (only `weight` and
  # `name` are handled by default, everything else is "ignore"), which is exactly
  # the silent-loss case this helper exists to prevent.
  graph <- igraph::as_undirected(
    graph, mode = "collapse",
    edge.attr.comb = list(.netkit_s = "sum", .netkit_d = "min", "ignore")
  )

  w$strength <- igraph::edge_attr(graph, ".netkit_s")
  w$distance <- igraph::edge_attr(graph, ".netkit_d")

  graph <- igraph::delete_edge_attr(graph, ".netkit_s")
  graph <- igraph::delete_edge_attr(graph, ".netkit_d")

  list(graph = graph, w = w)
}
