#' Calculate Network Roles Based on Within-Module Z-Score and Participation Coefficient
#'
#' Implements the node role classification system of Guimerà & Amaral (2005) by calculating
#' the within-module degree z-score and the participation coefficient for each node in a network.
#' Nodes are assigned to one of seven role categories (R1–R7) based on their local modular connectivity.
#'
#' If no community structure is provided, modules are automatically detected using the specified clustering method.
#' The function can optionally produce a 2D role plot (z vs. P) highlighting the canonical role regions.
#'
#' @param graph An `igraph` object representing the network.
#' @param communities Optional. A community clustering object (as returned by an `igraph` clustering function), or a named membership vector. If `NULL`, community detection is performed using `cluster.method`.
#' @param cluster.method Character. Clustering algorithm to use if `communities` is `NULL`. Default is `"spinglass"`. Passed to `find_modules()`.
#' @param plot Logical. Whether to generate a 2D plot of participation coefficient (P) vs. within-module z-score (z). Default is `TRUE`.
#' @param highlight_roles Logical. If `TRUE`, the role regions in the z–P plane are shaded for visual clarity. Default is `TRUE`.
#' @param hub_z Numeric. Threshold for defining hubs in terms of within-module z-score. Default is `2.5`.
#' @param label_region Optional character vector of role labels (e.g., `c("R4", "R7")`) indicating which role regions should have their nodes labeled in the plot. Default is `NULL`.
#' @param label.size Numeric. Base font size for plot text. Default is `12`.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{`plot`}{The `ggplot2` object (only if `plot = TRUE`).}
#'   \item{`roles_definitions`}{A data frame describing the seven role types and their conditions.}
#'   \item{`result`}{A data frame with node-level information: node name, module, z-score, participation coefficient, and assigned role.}
#' }
#'
#' @details
#' The node roles are defined as:
#'
#' \tabular{ll}{
#' R1 \tab Ultra-peripheral (non-hub): \eqn{z < 2.5, P <= 0.05} \cr
#' R2 \tab Peripheral (non-hub): \eqn{z < 2.5, 0.05 < P <= 0.6} \cr
#' R3 \tab Non-hub connector: \eqn{z < 2.5, 0.6 < P <= 0.8} \cr
#' R4 \tab Non-hub kinless: \eqn{z < 2.5, P > 0.8} \cr
#' R5 \tab Provincial hub: \eqn{z >= 2.5, P <= 0.3} \cr
#' R6 \tab Connector hub: \eqn{z >= 2.5, 0.3 < P <= 0.75} \cr
#' R7 \tab Kinless hub: \eqn{z >= 2.5, P > 0.75} \cr
#' }
#'
#' @references
#' Guimerà, R., & Amaral, L. A. N. (2005). Functional cartography of complex metabolic networks. *Nature*, 433(7028), 895–900. \doi{10.1038/nature03288}
#'
#' @seealso [find_modules()], [igraph::cluster_spinglass()], [igraph::membership()]
#'
#' @examples
#' \dontrun{
#' g <- igraph::sample_gnp(200, 0.05, directed = F)
#' igraph::V(g)$name <- as.character(1:200)
#' result <- calculate_roles(g, plot = TRUE)
#' head(result$result)
#' }
#'
#' @importFrom ggplot2 theme_bw ggplot annotate geom_point scale_x_continuous scale_color_manual labs
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
                            label.size = 12) {

  # Extract membership vector
  if (is.null(communities)) {

    modules <- find_modules(graph, method = cluster.method, plot = FALSE, return_subgraphs = FALSE)
    membership <- setNames(modules$module_table$module, modules$module_table$node)

  } else if ("communities" %in% class(communities)) {
    membership <- igraph::membership(communities)
  } else if (is.vector(communities) && length(communities) == igraph::vcount(graph)) {
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
                            paste0("z < ", hub_z, " & P <= 0.05"),
                            paste0("z < ", hub_z, " & 0.05 < P & P <= 0.6"),
                            paste0("z < ", hub_z, " & 0.6 < P & P <= 0.8"),
                            paste0("z < ", hub_z, " & P > 0.8"),
                            paste0("z >= ", hub_z, " & P <= 0.25"),
                            paste0("z >= ", hub_z, " & 0.25 < P & P <= 0.75"),
                            paste0("z >= ", hub_z, " & P > 0.75")
                          ))


  vnames <- names(membership)
  if (is.null(vnames)) {
    vnames <- as.character(seq_len(length(membership)))
  }

  degrees <- degree(graph, mode = "all")[vnames]

  # Initialize roles dataframe
  roles_df <- tibble::tibble(
    node = vnames,
    module = as.integer(membership[vnames]),
    z = NA_real_,
    p = NA_real_,
    role = NA_character_,
    stringsAsFactors = FALSE
  )

  # Compute within-module z-score
  for (mod in unique(roles_df$module)) {
    mod_nodes <- roles_df$node[roles_df$module == mod]
    if (length(mod_nodes) <= 1) {
      roles_df$z[roles_df$node %in% mod_nodes] <- 0
      next
    }

    subg <- induced_subgraph(graph, mod_nodes)
    ki <- degree(subg)
    mean_ki <- mean(ki)
    sd_ki <- sd(ki)

    z_vals <- if (is.na(sd_ki) || sd_ki == 0) rep(0, length(ki)) else (ki - mean_ki) / sd_ki
    roles_df$z[match(names(ki), roles_df$node)] <- z_vals
  }

  # Compute participation coefficient
  for (i in seq_len(nrow(roles_df))) {
    node <- roles_df$node[i]
    neighbors <- neighbors(graph, node, mode = "all")
    if (length(neighbors) == 0) {
      roles_df$p[i] <- 0
      next
    }

    neighbor_names <- V(graph)$name[neighbors]
    k_i <- degrees[node]
    neighbor_modules <- membership[neighbor_names]
    k_i_m <- table(neighbor_modules)
    sum_frac_sq <- sum((k_i_m / k_i)^2)
    roles_df$p[i] <- 1 - sum_frac_sq
  }

  # Classify roles based on z and p
  for (i in seq_len(nrow(roles_df))) {
    z <- roles_df$z[i]
    p <- roles_df$p[i]

    if (is.na(z) || is.na(p)) {
      roles_df$role[i] <- NA
    } else if (z < hub_z) {
      if (p <= 0.05) roles_df$role[i] <- "R1"
      else if (p <= 0.60) roles_df$role[i] <- "R2"
      else if (p <= 0.80) roles_df$role[i] <- "R3"
      else roles_df$role[i] <- "R4"
    } else {
      if (p <= 0.30) roles_df$role[i] <- "R5"
      else if (p <= 0.75) roles_df$role[i] <- "R6"
      else roles_df$role[i] <- "R7"
    }
  }


  # Plot if requested
  if (plot) {
    p <- ggplot(roles_df, aes(x = p, y = z, color = role))

    if (highlight_roles) {

      p <- p + annotate("rect",
                          xmin = -Inf, xmax = 0.05,
                          ymin = -Inf, ymax = hub_z,
                          alpha = 0.2, fill = "black") +
        annotate("rect",
                 xmin = 0.05, xmax = 0.6,
                 ymin = -Inf, ymax = hub_z,
                 alpha = 0.2, fill = "red") +
        annotate("rect",
                 xmin = 0.6, xmax = 0.8,
                 ymin = -Inf, ymax = hub_z,
                 alpha = 0.2, fill = "green") +
        annotate("rect",
                 xmin = 0.8, xmax = Inf,
                 ymin = -Inf, ymax = hub_z,
                 alpha = 0.2, fill = "darkblue") +
        annotate("rect",
                 xmin = -Inf, xmax = 0.25,
                 ymin = hub_z, ymax = Inf,
                 alpha = 0.2, fill = "yellow") +
        annotate("rect",
                 xmin = 0.25, xmax = 0.75,
                 ymin = hub_z, ymax = Inf,
                 alpha = 0.2, fill = "brown") +
        annotate("rect",
                 xmin = 0.75, xmax = Inf,
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

    return(list(
      plot = p,
      roles_definitions = roles_def,
      result = roles_df
    ))

  } else {

    return(list(
      roles_definitions = roles_def,
      result = roles_df
    ))
  }
}
