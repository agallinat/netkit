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
#'     `role` are `NA` for any vertex whose module, or whose neighbours' modules, are
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
#' memberships of its neighbours, so dropping a neighbour's module silently distorts
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
                            thresholds = NULL) {

  # Results are keyed by vertex name throughout, so names must exist.
  graph <- as_netkit_graph(graph, backfill_names = TRUE)

  # Single source of truth for the role boundaries; see R/utils-roles.R.
  th <- as_role_thresholds(thresholds)

  # Extract membership vector
  if (is.null(communities)) {

    # min_size = 1 is deliberate. find_modules() defaults to min_size = 3, which
    # discards small modules -- and a node whose module was discarded is absent
    # from module_table, so its membership is NA. That silently corrupted the
    # participation coefficient of *retained* neighbours, because the NA was
    # dropped from the neighbour tally while the full degree was still used as the
    # denominator. Roles are a per-node measure, so there is no reason to filter
    # modules by size here at all.
    modules <- find_modules(graph, method = cluster.method, min_size = 1,
                            plot = FALSE, return_subgraphs = FALSE)
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

  # Compute within-module z-score. NA modules are skipped rather than treated as a
  # module of their own: a membership vector containing NAs would otherwise select
  # NA node names and make induced_subgraph() fail with "Invalid vertex names".
  # Those nodes keep z = NA and so are classified with role = NA.
  for (mod in unique(roles_df$module[!is.na(roles_df$module)])) {
    idx <- which(roles_df$module == mod)
    mod_nodes <- roles_df$node[idx]
    if (length(mod_nodes) <= 1) {
      roles_df$z[idx] <- 0
      next
    }

    subg <- igraph::induced_subgraph(graph, mod_nodes)
    ki <- igraph::degree(subg)
    mean_ki <- mean(ki)
    sd_ki <- stats::sd(ki)

    z_vals <- if (is.na(sd_ki) || sd_ki == 0) rep(0, length(ki)) else (ki - mean_ki) / sd_ki
    roles_df$z[match(names(ki), roles_df$node)] <- z_vals
  }

  # Compute participation coefficient
  for (i in seq_len(nrow(roles_df))) {
    node <- roles_df$node[i]
    nbrs <- igraph::neighbors(graph, node, mode = "all")
    if (length(nbrs) == 0) {
      roles_df$p[i] <- 0
      next
    }

    neighbor_names <- igraph::V(graph)$name[nbrs]
    neighbor_modules <- membership[neighbor_names]
    k_i_m <- table(neighbor_modules)

    # The denominator must be the number of neighbours actually tallied, not the
    # node's full degree: table() drops neighbours with NA membership, so using
    # the full degree makes the fractions sum to less than 1 and inflates P.
    # These agree whenever every neighbour has a module, and stay well defined
    # when some do not.
    k_i <- sum(k_i_m)

    roles_df$p[i] <- if (k_i == 0) NA_real_ else 1 - sum((k_i_m / k_i)^2)
  }

  # Classify roles based on z and p
  for (i in seq_len(nrow(roles_df))) {
    z <- roles_df$z[i]
    p <- roles_df$p[i]

    if (is.na(z) || is.na(p)) {
      roles_df$role[i] <- NA
    } else if (z < hub_z) {
      if (p <= th[["R1_R2"]]) roles_df$role[i] <- "R1"
      else if (p <= th[["R2_R3"]]) roles_df$role[i] <- "R2"
      else if (p <= th[["R3_R4"]]) roles_df$role[i] <- "R3"
      else roles_df$role[i] <- "R4"
    } else {
      if (p <= th[["R5_R6"]]) roles_df$role[i] <- "R5"
      else if (p <= th[["R6_R7"]]) roles_df$role[i] <- "R6"
      else roles_df$role[i] <- "R7"
    }
  }


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
  return(list(
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
      "; hub z-score threshold = ", hub_z,
      "; participation boundaries R1/R2 = ", th[["R1_R2"]],
      ", R2/R3 = ", th[["R2_R3"]],
      ", R3/R4 = ", th[["R3_R4"]],
      ", R5/R6 = ", th[["R5_R6"]],
      ", R6/R7 = ", th[["R6_R7"]]
    ),
    roles_definitions = roles_def
  ))
}
