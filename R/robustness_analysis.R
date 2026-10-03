#' Network Robustness Analysis via Node Removal Simulation
#'
#' Simulates the removal of nodes from a network using various strategies and
#' evaluates how the structure degrades using selected robustness metrics.
#' Useful for assessing the vulnerability or resilience of a graph.
#'
#' @param graph An \code{igraph} object representing the network to plot or a
#'   data frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes. Must be undirected; directed graphs
#'   will be converted.
#' @param removal_strategy Character. Strategy used for node removal. Options are:
#'   \code{"random"}, \code{"degree"}, \code{"betweenness"}, or the name of a numeric vertex attribute.
#'   Custom attributes are interpreted as priority scores (higher = removed first).
#' @param steps Integer. Number of removal steps (default: 50).
#' @param metrics Character vector. Structural metrics to compute at each step.
#'   Options include: \code{"lcc_size"}, \code{"efficiency"}, and \code{"n_components"}.
#' @param n_reps Integer. Number of simulation repetitions (only relevant if \code{removal_strategy = "random"}).
#' @param plot Logical. If \code{TRUE}, a robustness plot is generated.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{plot}}{A \code{ggplot2} object showing the evolution of the selected
#'     metrics as nodes are progressively removed, or \code{NULL} when
#'     \code{plot = FALSE}. The element is always present, so the return shape does
#'     not depend on the arguments.}
#'   \item{\code{all_results}}{A data frame with simulation results across all steps and repetitions.}
#'   \item{\code{summary}}{A summarized data frame (mean and SD) if \code{n_reps > 1}, otherwise raw results.}
#'   \item{\code{auc}}{Named list of AUC (area under the curve) values for each selected metric.}
#' }
#'
#' @details
#' This function builds on classic approaches in network science for evaluating structural robustness,
#' simulating progressive node removal and quantifying the degradation of key topological features.
#'
#' For deterministic strategies (\code{"degree"}, \code{"betweenness"}, or custom attributes),
#' nodes are removed in a fixed priority order. For the \code{"random"} strategy, the process is
#' repeated \code{n_reps} times, and the results are aggregated.
#'
#' The available metrics are:
#' \itemize{
#'   \item \strong{Largest Connected Component}: size of the largest remaining component (\code{lcc_size}).
#'   \item \strong{Global Efficiency}: average inverse shortest path length among all pairs (\code{efficiency}).
#'   \item \strong{Number of Components}: total number of disconnected components (\code{n_components}).
#' }
#'
#' Additionally, Area Under the Curve (AUC) is calculated for each metric, providing a scalar summary
#' of robustness. A higher AUC indicates greater resilience (i.e., slower degradation).
#'
#' @references
#' Albert R, Jeong H, Barabási AL. Error and attack tolerance of complex networks.
#' Nature. 2000;406(6794):378–382. \doi{10.1038/35019019}
#'
#' @examples
#' g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
#'
#' # `steps` is reduced from its default of 50 to keep the example fast.
#' res <- robustness_analysis(g, removal_strategy = "degree", steps = 10,
#'                            metrics = c("lcc_size", "n_components"),
#'                            plot = FALSE, seed = 1)
#' res$auc
#' head(res$summary)
#'
#' @importFrom igraph is_igraph is_directed as_undirected vertex_attr_names V degree betweenness delete_vertices components vcount distances vertex_attr vertex_attr<-
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom dplyr bind_rows group_by summarise across all_of %>%
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal scale_color_manual
#' @importFrom pracma trapz
#' @importFrom stats sd
#'
#' @export
robustness_analysis <- function(graph,
                                removal_strategy = c("random", "degree", "betweenness"),
                                steps = 50,
                                metrics = c("lcc_size", "efficiency", "n_components"),
                                n_reps = 50,
                                plot = TRUE,
                                seed = NULL) {

  # Validate input
  graph <- as_netkit_graph(graph)

  if (is_directed(graph)) graph <- as_undirected(graph, mode = "collapse")

  # Resolve the unevaluated default (a length-3 vector) to its first option
  if (length(removal_strategy) > 1) {
    removal_strategy <- match.arg(removal_strategy, c("random", "degree", "betweenness"))
  }

  # Handle flexible strategy
  if (length(removal_strategy) == 1 && removal_strategy %in% c("random", "degree", "betweenness")) {
    strategy <- removal_strategy
    custom_vector <- NULL
  } else if (length(removal_strategy) == 1 && removal_strategy %in% vertex_attr_names(graph)) {
    strategy <- "custom"
    custom_vector <- vertex_attr(graph, removal_strategy)
  } else {
    stop("Invalid `removal_strategy`. Use 'random', 'degree', 'betweenness', or a valid vertex attribute name.")
  }

  metrics <- match.arg(metrics, several.ok = TRUE)
  set.seed(seed)

  n <- vcount(graph) - 1
  step_size <- ceiling(n / steps)
  result_list <- list()
  reps <- if (strategy == "random") n_reps else 1

  if (strategy == "random" && n_reps > 1) {
    pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
  }

  for (rep in 1:reps) {
    removal_order <- switch(strategy,
                            random = sample(V(graph)),
                            degree = V(graph)[order(-degree(graph))],
                            betweenness = V(graph)[order(-betweenness(graph))],
                            custom = V(graph)[order(-custom_vector)]
    )

    results <- list()

    for (i in seq(0, n, by = step_size)) {
      to_remove <- removal_order[1:i]
      g_tmp <- delete_vertices(graph, to_remove)

      row <- list(
        rep = rep,
        removed = length(to_remove),
        removed_frac = length(to_remove) / n
      )

      if ("lcc_size" %in% metrics) {
        comps <- components(g_tmp)
        row$lcc_size <- max(comps$csize)
      }

      if ("efficiency" %in% metrics) {
        n_tmp <- vcount(g_tmp)
        if (n_tmp < 2) {
          # Global efficiency averages 1/d over ordered pairs, and a graph with
          # fewer than two vertices has none: the n(n-1) denominator is 0, so
          # computing it directly yields NaN. With the default steps = 50 the
          # removal loop reaches a single remaining vertex for any graph of about
          # 51 nodes or fewer, and that NaN then propagated into the AUC. A
          # collapsed network has no paths left, so report 0.
          row$efficiency <- 0
        } else {
          sp <- distances(g_tmp)
          inv_sp <- 1 / sp
          inv_sp[is.infinite(inv_sp)] <- 0
          row$efficiency <- sum(inv_sp) / (n_tmp^2 - n_tmp)
        }
      }

      if ("n_components" %in% metrics) {
        comps <- components(g_tmp)
        row$n_components <- comps$no
      }

      results[[length(results) + 1]] <- row
    }

    result_list[[rep]] <- bind_rows(results)

    # update progress bar
    if (strategy == "random" && n_reps > 1) {
      setTxtProgressBar(pb, rep)
    }

  }

  if (strategy == "random" && n_reps > 1) {
    close(pb)
  }

  all_results <- bind_rows(result_list)

  # --- Aggregate if needed ---
  if (removal_strategy == "random" && n_reps > 1) {
    metric_cols <- metrics  # Only include selected metrics
    summary <- all_results %>%
      group_by(removed_frac) %>%
      summarise(
        across(all_of(metric_cols), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
        .groups = "drop"
      )
  } else {
    summary <- all_results
  }

  # --- AUC calculation ---
  auc_list <- list()
  for (metric in metrics) {
    metric_vec <- if (removal_strategy == "random" && n_reps > 1) {
      summary[[paste0(metric, "_mean")]]
    } else {
      summary[[metric]]
    }

    metric_vec <- normalize_metric(metric_vec)
    auc_val <- pracma::trapz(summary$removed_frac, metric_vec)
    auc_list[[metric]] <- auc_val
  }

  # --- Plot ---
  if (plot) {
    p <- ggplot(summary, aes(x = removed_frac)) +
      labs(x = "Fraction of Nodes Removed", y = NULL,
           title = paste("Network Robustness:", removal_strategy)) +
      theme_minimal(base_size = 13)

    if ("lcc_size" %in% metrics) {
      if ("lcc_size_mean" %in% names(summary)) {
        p <- p +
          geom_line(aes(y = normalize_metric(lcc_size_mean), color = "LCC Size"))
      } else {
        p <- p + geom_line(aes(y = normalize_metric(lcc_size), color = "LCC Size"))
      }
    }

    if ("efficiency" %in% metrics) {
      if ("efficiency_mean" %in% names(summary)) {
        p <- p +
          geom_line(aes(y = normalize_metric(efficiency_mean), color = "Efficiency"))
      } else {
        p <- p + geom_line(aes(y = normalize_metric(efficiency), color = "Efficiency"))
      }
    }

    if ("n_components" %in% metrics) {
      if ("n_components_mean" %in% names(summary)) {
        p <- p +
          geom_line(aes(y = normalize_metric(n_components_mean), color = "Components"))
      } else {
        p <- p + geom_line(aes(y = normalize_metric(n_components), color = "Components"))
      }
    }

    p <- p +
      scale_color_manual(values = c("LCC Size" = "steelblue", "Efficiency" = "darkgreen", "Components" = "red")) +
      # No fill aesthetic is mapped, so naming one made ggplot2 report
      # "Ignoring unknown labels: fill".
      labs(color = "Metric")

  } else {

    p <- NULL

  }

  # `plot` is always present, and NULL when plot = FALSE, so that the return
  # shape does not depend on the arguments.
  return(list(
    plot = p,
    all_results = all_results,
    summary = summary,
    auc = auc_list
  ))
}
