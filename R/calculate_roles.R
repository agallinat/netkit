#' Calculate Network Roles Based on Within-Module Z-Score and Participation Coefficient
#'
#' Implements the node role classification system of Guimerà & Amaral (2005) by calculating
#' the within-module degree z-score and the participation coefficient for each node in a network.
#' Nodes are assigned to one of seven role categories (R1–R7) based on their local modular connectivity.
#'
#' If no community structure is provided, modules are automatically detected using the specified clustering method.
#' The function can optionally produce a 2D role plot (z vs. P) highlighting the canonical role regions.
#'
#' @param graph An `igraph` object representing the network, or a data frame
#'   containing a symbolic edge list in the first two columns. Additional columns are
#'   considered as edge attributes.
#' @param communities Optional. A community clustering object (as returned by an `igraph` clustering function), or a named membership vector. If `NULL`, community detection is performed using `cluster.method`.
#' @param cluster.method Character. Clustering algorithm to use if `communities` is `NULL`. Default is `"spinglass"`. Passed to `find_modules()`.
#' @param plot Logical. Whether to generate a 2D plot of participation coefficient (P) vs. within-module z-score (z). Default is `TRUE`.
#' @param highlight_roles Logical. If `TRUE`, the role regions in the z–P plane are shaded for visual clarity. Default is `TRUE`.
#' @param hub_z Numeric. Threshold for defining hubs in terms of within-module z-score. Default is `2.5`.
#' @param label_region Optional character vector of role labels (e.g., `c("R4", "R7")`) indicating which role regions should have their nodes labeled in the plot. Default is `NULL`.
#' @param label.size Numeric. Base font size for plot text. Default is `12`.
#' @param thresholds Optional named numeric vector overriding one or more of the
#'   participation-coefficient boundaries between roles. Names must be drawn from
#'   `R1_R2`, `R2_R3`, `R3_R4`, `R5_R6` and `R6_R7`; unnamed entries, unknown
#'   names, values outside `[0, 1]` and non-monotonic sets are rejected. `NULL`
#'   (default) uses the published values -- see Details.
#'
#' @return A list with five elements:
#' \describe{
#'   \item{`plot`}{A `ggplot2` object, or `NULL` when `plot = FALSE`. The element is
#'     always present, so the return shape does not depend on the arguments.}
#'   \item{`result`}{A data frame with node-level information: node name, module,
#'     z-score, participation coefficient, and assigned role. It has one row per
#'     graph vertex. The exception is `cluster.method = "spinglass"`, which can only
#'     be run on the largest connected component of a disconnected graph; vertices
#'     outside it have no module and are absent, which raises a warning. `z`, `p` and
#'     `role` are `NA` for any vertex whose module, or whose neighbors' modules, are
#'     unknown.}
#'   \item{`graph`}{The input graph with `module`, `role_z`, `role_p` and `role`
#'     attached as vertex attributes, so that the classification can be passed
#'     straight to [plot_Net()] or [robustness_analysis()].}
#'   \item{`method`}{A human-readable description of the community detection used
#'     and the thresholds actually applied.}
#'   \item{`roles_definitions`}{A data frame describing the seven role types and
#'     their conditions, generated from the same thresholds the classifier used.}
#' }
#'
#' @details
#' When `communities` is `NULL`, community detection is delegated to
#' [find_modules()] with `min_size = 1`, so no module is discarded for being small
#' and every vertex receives a role. This matters for correctness as well as
#' coverage: the participation coefficient of a node is computed from the module
#' memberships of its neighbors, so dropping a neighbor's module silently distorts
#' the coefficient of the node that remains.
#'
#' The node roles are defined as follows, where `hub_z` defaults to 2.5 and the
#' participation-coefficient boundaries are those published in Guimera & Amaral
#' (2005):
#'
#' \tabular{ll}{
#' R1 \tab Ultra-peripheral (non-hub): \eqn{z < 2.5, P <= 0.05} \cr
#' R2 \tab Peripheral (non-hub): \eqn{z < 2.5, 0.05 < P <= 0.62} \cr
#' R3 \tab Non-hub connector: \eqn{z < 2.5, 0.62 < P <= 0.80} \cr
#' R4 \tab Non-hub kinless: \eqn{z < 2.5, P > 0.80} \cr
#' R5 \tab Provincial hub: \eqn{z >= 2.5, P <= 0.30} \cr
#' R6 \tab Connector hub: \eqn{z >= 2.5, 0.30 < P <= 0.75} \cr
#' R7 \tab Kinless hub: \eqn{z >= 2.5, P > 0.75} \cr
#' }
#'
#' Those five numbers are held in one place internally and are used by the
#' classifier, by the `roles_definitions` table and by the shaded bands of the
#' diagnostic plot alike, so the three cannot disagree. Override them with
#' `thresholds` if a different convention is wanted.
#'
#' @references
#' Guimerà, R., & Amaral, L. A. N. (2005). Functional cartography of complex metabolic networks. *Nature*, 433(7028), 895–900. \doi{10.1038/nature03288}
#'
#' @seealso [find_modules()], [igraph::cluster_spinglass()], [igraph::membership()]
#'
#' @examples
#' g <- igraph::sample_gnp(80, 0.08, directed = FALSE)
#' igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
#'
#' # "louvain" is used here rather than the "spinglass" default. spinglass cannot
#' # run on a disconnected graph, so find_modules() falls back to the largest
#' # connected component and the remaining nodes receive no role.
#' result <- calculate_roles(g, cluster.method = "louvain", plot = FALSE)
#' head(result$result)
#' result$roles_definitions
#'
#' @param weights Optional edge weights: `NULL` (default) to ignore them, the
#'   name of an edge attribute, or a numeric vector of length
#'   `igraph::ecount(graph)`. When supplied, the within-module z-score is
#'   computed from strengths and the participation coefficient from summed edge
#'   strengths per module, which is the weighted generalization given in
#'   Guimera & Amaral's supplementary material. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`. See
#'   [netkit-weights].
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @importFrom ggplot2 theme_bw ggplot aes annotate geom_point scale_x_continuous scale_color_manual labs theme
#' @importFrom ggrepel geom_text_repel
#'
#' @export
calculate_roles <- function(graph,
                            communities = NULL,
                            cluster.method = "spinglass",
                            plot = TRUE,
                            highlight_roles = TRUE,
                            hub_z = 2.5,
                            label_region = NULL,
                            label.size = 12,
                            thresholds = NULL,
                            weights = NULL,
                            weight_type = c("strength", "distance")) {

  # Results are keyed by vertex name throughout, so names must exist.
  graph <- as_netkit_graph(graph, backfill_names = TRUE)
  w <- as_netkit_weights(graph, weights, weight_type)

  # Single source of truth for the role boundaries; see R/utils-roles.R.
  th <- as_role_thresholds(thresholds)

  # Extract membership vector
  if (is.null(communities)) {

    # min_size = 1 is deliberate. find_modules() defaults to min_size = 3, which
    # discards small modules -- and a node whose module was discarded is absent
    # from module_table, so its membership is NA. That silently corrupted the
    # participation coefficient of *retained* neighbors, because the NA was
    # dropped from the neighbor tally while the full degree was still used as the
    # denominator. Roles are a per-node measure, so there is no reason to filter
    # modules by size here at all.
    modules <- find_modules(graph, method = cluster.method, min_size = 1,
                            plot = FALSE, return_subgraphs = FALSE,
                            weights = weights, weight_type = weight_type)
    membership <- stats::setNames(modules$module_table$module, modules$module_table$node)

  } else if ("communities" %in% class(communities)) {
    membership <- igraph::membership(communities)
  } else if (is.atomic(communities) && length(communities) == igraph::vcount(graph)) {
    # is.atomic() rather than is.vector(): igraph::membership() returns a named
    # vector carrying a "membership" class attribute, which is.vector() rejects.
    membership <- communities
  } else {
    stop("Communities must be an igraph clustering object, a membership vector or NULL to find modules.")
  }

  # Roles definition
  roles_def <- data.frame(Name = c("R1","R2","R3","R4","R5","R6","R7"),
                          Description = c("Ultra-peripheral (non-hub)",
                                          "Peripheral (non-hub)",
                                          "Non-hub connector",
                                          "Non-hub kinless",
                                          "Provincial hub",
                                          "Connector hub",
                                          "Kinless hub"),
                          Condition = c(
                            paste0("z < ", hub_z, " & P <= ", th[["R1_R2"]]),
                            paste0("z < ", hub_z, " & ", th[["R1_R2"]], " < P & P <= ", th[["R2_R3"]]),
                            paste0("z < ", hub_z, " & ", th[["R2_R3"]], " < P & P <= ", th[["R3_R4"]]),
                            paste0("z < ", hub_z, " & P > ", th[["R3_R4"]]),
                            paste0("z >= ", hub_z, " & P <= ", th[["R5_R6"]]),
                            paste0("z >= ", hub_z, " & ", th[["R5_R6"]], " < P & P <= ", th[["R6_R7"]]),
                            paste0("z >= ", hub_z, " & P > ", th[["R6_R7"]])
                          ))


  # An unnamed membership vector is positional, in vertex order, which is what
  # igraph::membership() returns. Naming it from the graph is required: inventing
  # index names ("1", "2", ...) instead meant every later vertex lookup failed with
  # "Invalid vertex names" whenever the graph's own names were anything else.
  if (is.null(names(membership))) {
    names(membership) <- igraph::V(graph)$name
  }
  vnames <- names(membership)

  # Warn rather than silently return a partial table. With communities = NULL this
  # should now cover every vertex; the exception is cluster.method = "spinglass",
  # which find_modules() can only run on the largest connected component.
  missing_nodes <- setdiff(igraph::V(graph)$name, vnames)
  if (length(missing_nodes) > 0) {
    warning(sprintf(
      "%d of %d vertices have no module assignment and are absent from the result: %s%s",
      length(missing_nodes), igraph::vcount(graph),
      paste(utils::head(missing_nodes, 5), collapse = ", "),
      if (length(missing_nodes) > 5) ", ..." else ""
    ))
  }

  # Initialize roles dataframe
  roles_df <- tibble::tibble(
    node = vnames,
    module = as.integer(membership[vnames]),
    z = NA_real_,
    p = NA_real_,
    role = NA_character_
  )

  # Compute within-module z-score.
  #
  # Modules indexed by vertex position, as small integer codes, so every tally
  # below is an arithmetic pass rather than a lookup by name. NA membership stays
  # NA: such vertices keep z = NA and so are classified with role = NA. Treating
  # them as a module of their own would be wrong, and the previous code could not
  # have done so anyway -- it selected node names by module and an NA name makes
  # induced_subgraph() fail with "Invalid vertex names".
  gnames <- igraph::V(graph)$name
  n_v <- igraph::vcount(graph)

  mod_by_pos <- membership[match(gnames, names(membership))]
  mod_levels <- sort(unique(mod_by_pos[!is.na(mod_by_pos)]))
  mod_code <- match(mod_by_pos, mod_levels)
  n_mod <- length(mod_levels)

  # Within-module degree in one pass over the edge list, rather than one
  # induced_subgraph() per module. The subgraph approach is O(modules x (V + E)):
  # with min_size = 1 a 6000-vertex graph can carry thousands of modules, and
  # induced_subgraph() then accounted for the bulk of the runtime -- 1.32 s at
  # 3000 modules against 0.07 s at 10. An edge contributes to its endpoints'
  # within-module degree exactly when both ends share a module, which is the same
  # thing induced_subgraph() computes by building the subgraph first.
  #
  # Self-loops and parallel edges land where igraph::degree(loops = TRUE) and
  # igraph::strength() put them: a self-loop adds its contribution to the vertex
  # twice, once per endpoint, and parallel edges count separately.
  el <- igraph::as_edgelist(graph, names = FALSE)
  e_from <- el[, 1L]
  e_to <- el[, 2L]
  # Both endpoints need the non-NA guard, not just one: `TRUE & NA` is NA, so
  # testing only `e_from` let an edge to an unknown-module neighbor put an NA
  # into `internal` and from there into a subscript.
  internal <- !is.na(mod_code[e_from]) & !is.na(mod_code[e_to]) &
    mod_code[e_from] == mod_code[e_to]
  e_val <- if (w$weighted) w$strength else rep.int(1, igraph::ecount(graph))

  within_deg <- numeric(n_v)
  if (any(internal)) {
    contrib <- c(e_val[internal], e_val[internal])
    ends_v <- c(e_from[internal], e_to[internal])
    acc <- rowsum(contrib, ends_v, reorder = FALSE)
    within_deg[as.integer(rownames(acc))] <- acc[, 1L]
  }

  # Per-module mean and sample sd, computed in two passes so the variance is a
  # sum of squared deviations rather than E[x^2] - E[x]^2. The one-pass form
  # cancels badly, and z feeds a hub threshold, so the difference would be a
  # misclassification rather than a rounding artifact.
  valid <- !is.na(mod_code)
  code_v <- mod_code[valid]
  cnt <- tabulate(code_v, nbins = max(n_mod, 1L))

  sum_m <- numeric(max(n_mod, 1L))
  if (length(code_v)) {
    acc <- rowsum(within_deg[valid], code_v, reorder = FALSE)
    sum_m[as.integer(rownames(acc))] <- acc[, 1L]
  }
  mean_m <- ifelse(cnt > 0, sum_m / cnt, NA_real_)

  dev <- within_deg - mean_m[mod_code]
  ss_m <- numeric(max(n_mod, 1L))
  if (length(code_v)) {
    acc <- rowsum(dev[valid]^2, code_v, reorder = FALSE)
    ss_m[as.integer(rownames(acc))] <- acc[, 1L]
  }
  # cnt == 1 gives NA, which the zero-variance branch below turns into z = 0 --
  # the same answer the old `length(mod_nodes) <= 1` short-circuit gave.
  sd_m <- ifelse(cnt > 1, sqrt(ss_m / (cnt - 1)), NA_real_)

  sd_v <- sd_m[mod_code]
  z_by_pos <- ifelse(is.na(sd_v) | sd_v == 0, 0, dev / sd_v)
  z_by_pos[!valid] <- NA_real_

  roles_df$z <- unname(z_by_pos)[match(roles_df$node, gnames)]

  # Compute participation coefficient.
  #
  # One aggregation over the whole graph rather than one pass per vertex. The
  # per-vertex loop this replaces called neighbors(graph, node) with a character
  # `node`, so igraph resolved a name on every iteration, and then evaluated
  # `igraph::V(graph)$name[nbrs]`, which materializes the entire name vector to
  # take a handful of elements from it. Both are inside the loop, so the work grew
  # as the square of the vertex count: at 6000 vertices, vertex.attributes() alone
  # was 42% of the runtime and the name machinery about 80% of it between them.
  #
  # The neighbor and incident-edge lists come from adjacent_vertices() and
  # incident_edges(), which are the batched forms of the neighbors() and
  # incident() calls used before. Using the same primitives keeps the treatment of
  # self-loops and parallel edges identical rather than re-deriving it.
  #
  # `gnames`, `mod_code` and `n_mod` are the ones built for the z-score above:
  # modules as small integer codes indexed by vertex position, so a neighbor's
  # module is an integer lookup and can index a matrix column directly.
  if (w$weighted) {
    # Pair each incident edge with its *own* far endpoint. Zipping incident()
    # against neighbors() would not do: both are sorted, but by edge id and by
    # vertex id respectively, so a mis-pairing would attribute one neighbor's
    # strength to another neighbor's module and still produce a coefficient in
    # [0, 1].
    inc <- igraph::incident_edges(graph, igraph::V(graph), mode = "all")
    e_of_pair <- as.integer(unlist(inc, use.names = FALSE))
    owner <- rep.int(seq_len(n_v), lengths(inc))

    ea <- el[e_of_pair, 1L]
    eb <- el[e_of_pair, 2L]
    far <- ifelse(ea == owner, eb, ea)

    pair_mod <- mod_code[far]
    pair_val <- w$strength[e_of_pair]
    deg_all <- lengths(inc)
  } else {
    nbr <- igraph::adjacent_vertices(graph, igraph::V(graph), mode = "all")
    pair_mod <- mod_code[as.integer(unlist(nbr, use.names = FALSE))]
    owner <- rep.int(seq_len(n_v), lengths(nbr))
    pair_val <- rep.int(1, length(owner))
    deg_all <- lengths(nbr)
  }

  # Neighbors whose module is unknown are dropped, so the numerator and the
  # denominator are both derived from the surviving pairs -- the same place. The
  # denominator must not be the node's full degree: dropping a neighbor from the
  # tally while still dividing by the full degree makes the fractions sum to less
  # than 1 and inflates P, which is the defect this function shipped once before.
  keep <- !is.na(pair_mod)

  # Duplicate (vertex, module) entries are summed by sparseMatrix, which is what
  # makes this a tally. Unweighted that counts neighbors per module; weighted it
  # sums the strengths of the edges reaching each module -- the weighted
  # generalization, under which a node tied to one module by a strong edge and to
  # another by a weak one is not a connector, though counting alone would say so.
  k_im <- Matrix::sparseMatrix(
    i = owner[keep], j = pair_mod[keep], x = pair_val[keep],
    dims = c(n_v, max(n_mod, 1L))
  )

  k_i <- Matrix::rowSums(k_im)
  sum_sq <- Matrix::rowSums(k_im * k_im)

  # Zero-tally modules contribute nothing to the sum of squares, so restricting
  # to the modules a node actually touches is unnecessary: 1 - sum((k_im/k_i)^2)
  # expands to 1 - sum_sq / k_i^2.
  # unname(): lengths() carries the vertex names and ifelse() propagates the
  # attributes of its first argument, which would otherwise make `p` a named
  # column rather than a plain numeric one.
  p_by_pos <- unname(ifelse(deg_all == 0, 0,
                            ifelse(k_i == 0, NA_real_, 1 - sum_sq / k_i^2)))

  roles_df$p <- p_by_pos[match(roles_df$node, gnames)]

  # Classify roles based on z and p.
  #
  # Vectorized over the table rather than one pass per row. findInterval() encodes
  # the same ladder of comparisons: the boundaries are upper-inclusive, so
  # rightmost.closed is wrong here and the interval index is taken on `-p` against
  # the negated, reversed breaks -- which is what makes `p <= b` rather than
  # `p < b` the test, matching the published definition and the shaded plot bands.
  nonhub <- c(th[["R1_R2"]], th[["R2_R3"]], th[["R3_R4"]])
  hub <- c(th[["R5_R6"]], th[["R6_R7"]])

  pick <- function(p, breaks, labels) {
    labels[length(breaks) + 1L - findInterval(-p, -rev(breaks))]
  }

  is_hub <- !is.na(roles_df$z) & roles_df$z >= hub_z
  role <- rep(NA_character_, nrow(roles_df))
  ok <- !is.na(roles_df$z) & !is.na(roles_df$p)

  idx <- ok & !is_hub
  role[idx] <- pick(roles_df$p[idx], nonhub, c("R1", "R2", "R3", "R4"))
  idx <- ok & is_hub
  role[idx] <- pick(roles_df$p[idx], hub, c("R5", "R6", "R7"))
  roles_df$role <- role


  # Plot if requested
  if (plot) {
    p <- ggplot(roles_df, aes(x = p, y = z, color = role))

    if (highlight_roles) {

      # Band edges come from the same `th` the classifier uses, so a shaded region
      # can no longer disagree with the classification it illustrates.
      p <- p + annotate("rect",
                          xmin = -Inf, xmax = th[["R1_R2"]],
                          ymin = -Inf, ymax = hub_z,
                          alpha = 0.2, fill = "black") +
        annotate("rect",
                 xmin = th[["R1_R2"]], xmax = th[["R2_R3"]],
                 ymin = -Inf, ymax = hub_z,
                 alpha = 0.2, fill = "red") +
        annotate("rect",
                 xmin = th[["R2_R3"]], xmax = th[["R3_R4"]],
                 ymin = -Inf, ymax = hub_z,
                 alpha = 0.2, fill = "green") +
        annotate("rect",
                 xmin = th[["R3_R4"]], xmax = Inf,
                 ymin = -Inf, ymax = hub_z,
                 alpha = 0.2, fill = "darkblue") +
        annotate("rect",
                 xmin = -Inf, xmax = th[["R5_R6"]],
                 ymin = hub_z, ymax = Inf,
                 alpha = 0.2, fill = "yellow") +
        annotate("rect",
                 xmin = th[["R5_R6"]], xmax = th[["R6_R7"]],
                 ymin = hub_z, ymax = Inf,
                 alpha = 0.2, fill = "brown") +
        annotate("rect",
                 xmin = th[["R6_R7"]], xmax = Inf,
                 ymin = hub_z, ymax = Inf,
                 alpha = 0.2, fill = "gray")
    }

    p <- p +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(values = c("R1"="black",
                                    "R2"="red",
                                    "R3"="green",
                                    "R4"="darkblue",
                                    "R5"="orange",
                                    "R6"="brown",
                                    "R7"="darkgray")
                         ) +
      labs(x = "P", y = "z") +
      scale_x_continuous(breaks = seq.default(0, 1, 0.2),
                         limits = c(0,1)) +
      theme_bw(base_size = label.size) +
      theme(legend.position = "bottom")

    if (!is.null(label_region)) {

      p <- p +  geom_text_repel(data = roles_df[roles_df$role %in% label_region, ],
                      aes(label = node), color = "black",
                      point.padding =  3,
                      min.segment.length = 2,
                      size = 3)  # Names size,
    }

  } else {

    p <- NULL

  }

  # Annotate the graph so the classification chains into plot_Net(),
  # robustness_analysis() and the rest, as every other analysis function does.
  # match() keys by name and leaves NA for any vertex absent from roles_df, which
  # is the spinglass/largest-component case warned about above.
  vmap <- match(igraph::V(graph)$name, roles_df$node)
  igraph::vertex_attr(graph, "module") <- roles_df$module[vmap]
  igraph::vertex_attr(graph, "role_z") <- roles_df$z[vmap]
  igraph::vertex_attr(graph, "role_p") <- roles_df$p[vmap]
  igraph::vertex_attr(graph, "role")   <- roles_df$role[vmap]

  # `plot` is always present, and NULL when plot = FALSE, so that the return
  # shape does not depend on the arguments.
  return(as_netkit_result(list(
    plot = p,
    result = roles_df,
    graph = graph,
    method = paste0(
      "Guimera-Amaral roles from ",
      if (is.null(communities)) {
        paste0("modules detected by '", cluster.method, "'")
      } else {
        "a caller-supplied community structure"
      },
      " (", describe_weights(w), ")",
      "; hub z-score threshold = ", hub_z,
      "; participation boundaries R1/R2 = ", th[["R1_R2"]],
      ", R2/R3 = ", th[["R2_R3"]],
      ", R3/R4 = ", th[["R3_R4"]],
      ", R5/R6 = ", th[["R5_R6"]],
      ", R6/R7 = ", th[["R6_R7"]]
    ),
    roles_definitions = roles_def
  ), "calculate_roles"))
}
