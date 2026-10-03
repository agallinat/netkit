#' Generate Null-Model Graphs for Significance Testing
#'
#' Produces an ensemble of random graphs matched to an observed graph in some
#' respect, so that an observed metric can be compared against what that respect
#' alone would predict.
#'
#' This is what turns [summarize_graph_metrics()] from descriptive into
#' inferential. On its own, a modularity of 0.42 says nothing: random graphs with
#' the same degree sequence routinely reach 0.3. Pass an ensemble from this
#' function to [metric_significance()] to get the z-score and empirical p-value.
#'
#' @param graph An `igraph` object representing the network to analyze, or a data
#'   frame containing a symbolic edge list in the first two columns. Additional
#'   columns are considered as edge attributes.
#' @param model Which aspect of the observed graph to preserve:
#'   \describe{
#'     \item{`"rewire"`}{Degree-preserving edge rewiring (the default). Preserves
#'       the degree sequence exactly while staying close to the observed graph.}
#'     \item{`"configuration"`}{Draws a fresh graph with the same degree sequence
#'       via [igraph::sample_degseq()]. Preserves the degree sequence exactly but
#'       is otherwise unconstrained.}
#'     \item{`"erdos_renyi"`}{Matches only the vertex and edge counts. The
#'       weakest null, included as the deliberate contrast: comparing against it
#'       tells you whether an effect needs anything beyond density.}
#'   }
#' @param n Integer. Number of null graphs to generate. Default is `100`.
#' @param shuffle_weights Logical. If `TRUE` (default) and the observed graph
#'   carries a `weight` edge attribute, the observed weights are permuted across
#'   the null graph's edges, so the weight *distribution* is preserved while its
#'   placement is randomized. If `FALSE`, the null graphs carry no weights at
#'   all -- including under `model = "rewire"`, which would otherwise inherit the
#'   observed attribute bound to edges that have since been rewired.
#' @param seed Integer or `NULL`. Seed for reproducibility. When `NULL` (default)
#'   the caller's random stream is used and left alone.
#' @param niter Integer or `NULL`. Number of rewiring trials for
#'   `model = "rewire"`. When `NULL` (default), `10 * ecount(graph)`, which is the
#'   usual rule of thumb for mixing.
#'
#' @return A list of `n` `igraph` objects, with class `c("netkit_null", "list")`
#'   and the model recorded in its `"model"` attribute.
#'
#' @details
#' `"rewire"` and `"configuration"` both preserve the degree sequence, which is
#' almost always the right null for a topological claim: nearly every network
#' metric is partly determined by the degree sequence, so a null that does not
#' hold it fixed mostly rediscovers that the graph is heavy-tailed.
#'
#' They differ in how far they travel from the observed graph. `"rewire"` makes
#' local double-edge swaps, so the ensemble is centred on the observed graph and
#' is the more conservative choice. `"configuration"` resamples from scratch.
#'
#' `"configuration"` is attempted with `method = "vl"`, which produces simple
#' connected graphs, and falls back to `"configuration.simple"` with a warning
#' where that is not possible -- `"vl"` requires a connected realisation of the
#' degree sequence to exist, which fails for example when the graph has isolated
#' vertices.
#'
#' @references
#' Maslov, S., & Sneppen, K. (2002). Specificity and stability in topology of
#' protein networks. *Science*, 296(5569), 910-913.
#' \doi{10.1126/science.1065103}
#'
#' Newman, M. E. J., Strogatz, S. H., & Watts, D. J. (2001). Random graphs with
#' arbitrary degree distributions and their applications. *Physical Review E*,
#' 64(2), 026118. \doi{10.1103/PhysRevE.64.026118}
#'
#' @seealso [metric_significance()] to score an observed graph against the
#'   ensemble, [small_worldness()] for the classic application.
#'
#' @examples
#' g <- igraph::sample_pa(60, power = 1.5, directed = FALSE)
#'
#' # `n` is small here to keep the example fast; use 100 or more in practice.
#' nulls <- null_model(g, model = "rewire", n = 10, seed = 1)
#' length(nulls)
#'
#' # The degree sequence is preserved exactly, which is the point.
#' identical(sort(igraph::degree(nulls[[1]])), sort(igraph::degree(g)))
#'
#' @importFrom igraph rewire keeping_degseq sample_degseq sample_gnm degree
#' @importFrom igraph vcount ecount is_directed simplify as_undirected
#' @importFrom igraph edge_attr edge_attr<- edge_attr_names delete_edge_attr
#' @importFrom igraph vertex_attr vertex_attr<- V E as_edgelist edge_density
#'
#' @export
null_model <- function(graph,
                       model = c("rewire", "configuration", "erdos_renyi"),
                       n = 100,
                       shuffle_weights = TRUE,
                       seed = NULL,
                       niter = NULL) {

  model <- match.arg(model)

  graph <- as_netkit_graph(graph)

  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("'n' must be a single number of at least 1.", call. = FALSE)
  }
  n <- as.integer(n)

  # Only seed when asked: set.seed(NULL) would re-seed from the clock and
  # destroy the caller's stream, which is the defect robustness_analysis() had.
  if (!is.null(seed)) set.seed(seed)

  if (is_directed(graph)) {
    graph <- as_undirected(graph, mode = "collapse")
    message("Input graph converted to undirected for null-model generation.")
  }

  deg <- igraph::degree(graph)
  nv <- vcount(graph)
  ne <- ecount(graph)

  if (ne == 0) {
    stop("Cannot build a null model for a graph with no edges.", call. = FALSE)
  }

  if (is.null(niter)) niter <- 10 * ne

  weights <- if (shuffle_weights && "weight" %in% igraph::edge_attr_names(graph)) {
    igraph::edge_attr(graph, "weight")
  } else {
    NULL
  }

  # "vl" yields simple connected graphs but needs a connected realisation of the
  # degree sequence to exist; it fails outright on, for example, isolated
  # vertices. Probed once rather than inside the loop.
  degseq_method <- "vl"
  if (model == "configuration") {
    probe <- tryCatch(
      {
        igraph::sample_degseq(deg, method = "vl")
        TRUE
      },
      error = function(e) FALSE
    )
    if (!probe) {
      degseq_method <- "configuration.simple"
      warning("sample_degseq(method = \"vl\") cannot realise this degree ",
              "sequence as a connected simple graph (isolated vertices, or no ",
              "connected realisation). Falling back to ",
              "\"configuration.simple\", which may leave the null graphs ",
              "disconnected.", call. = FALSE)
    }
  }

  nulls <- vector("list", n)
  for (i in seq_len(n)) {
    null_g <- switch(
      model,
      rewire = igraph::rewire(graph, with = igraph::keeping_degseq(niter = niter)),
      configuration = igraph::sample_degseq(deg, method = degseq_method),
      erdos_renyi = igraph::sample_gnm(nv, ne, directed = FALSE)
    )

    # sample_degseq() and sample_gnm() return unnamed graphs. Names are needed by
    # anything that keys results by vertex, so carry the observed ones over --
    # the degree sequence is preserved positionally, so this is consistent.
    if (is.null(igraph::vertex_attr(null_g, "name")) &&
        !is.null(igraph::vertex_attr(graph, "name"))) {
      igraph::vertex_attr(null_g, "name") <- igraph::vertex_attr(graph, "name")
    }

    if (!is.null(weights)) {
      # Permute rather than copy: the distribution of weights is preserved, but
      # which edge carries which is randomized, so a weighted metric is tested
      # against weight *placement* and not merely against the weight values.
      igraph::edge_attr(null_g, "weight") <-
        sample(weights, ecount(null_g), replace = ecount(null_g) > length(weights))
    } else if ("weight" %in% igraph::edge_attr_names(null_g)) {
      # model = "rewire" starts from the observed graph, so it inherits the
      # observed `weight` attribute -- now bound to edges that have been rewired.
      # That is neither the observed weighting nor a declared permutation of it,
      # so leaving it in place would make shuffle_weights = FALSE mean something
      # undocumented. Drop it: the nulls are then unambiguously unweighted.
      null_g <- igraph::delete_edge_attr(null_g, "weight")
    }

    nulls[[i]] <- null_g
  }

  structure(nulls, class = c("netkit_null", "list"), model = model)
}

#' Score an Observed Graph Against a Null Ensemble
#'
#' Compares each of a graph's global metrics against the distribution that metric
#' takes over a null ensemble, and reports the z-score and empirical p-value.
#'
#' @param graph An `igraph` object or a data frame edge list, as elsewhere in
#'   netkit.
#' @param null A list of null graphs, as produced by [null_model()]. If `NULL`
#'   (default), an ensemble is generated internally with `n_null` graphs under
#'   `model`.
#' @param metrics Character vector of metric names to test, matching columns of
#'   [summarize_graph_metrics()]. If `NULL` (default), every numeric metric that
#'   varies across the ensemble is used.
#' @param n_null Integer. Number of null graphs when `null` is `NULL`. Default is
#'   `100`.
#' @param model Passed to [null_model()] when `null` is `NULL`.
#' @param weights,weight_type Passed to [summarize_graph_metrics()]. See
#'   [netkit-weights].
#' @param seed Integer or `NULL`. Seed for reproducibility.
#' @param plot Logical. Whether to build a diagnostic plot. Default is `TRUE`.
#' @param label.size Numeric. Base font size for plot text. Default is `12`.
#'
#' @return A list with:
#' \describe{
#'   \item{`plot`}{A faceted `ggplot2` object showing each metric's null
#'     distribution with the observed value marked, or `NULL` when
#'     `plot = FALSE`. The element is always present.}
#'   \item{`result`}{A tibble with one row per metric: `metric`, `observed`,
#'     `null_mean`, `null_sd`, `z`, `p_empirical`, `ci_lower` and `ci_upper`
#'     (the 2.5th and 97.5th percentiles of the null).}
#'   \item{`graph`}{The input graph, unchanged, so the call still chains.}
#'   \item{`method`}{A human-readable description of the test performed.}
#' }
#'
#' @details
#' `p_empirical` is two-sided and uses the `(r + 1) / (n + 1)` convention, where
#' `r` counts null values at least as extreme as the observed one. It is
#' therefore never exactly zero: with `n` nulls the smallest attainable p-value
#' is `1 / (n + 1)`, so testing against 20 nulls cannot produce evidence at
#' `p < 0.05` however large the effect. Choose `n_null` accordingly.
#'
#' `z` is `NA` when the null distribution has zero variance -- which happens
#' legitimately for metrics the null model holds fixed, such as `Nodes`, `Edges`
#' and `Density` under a degree-preserving model. Those metrics are excluded from
#' the automatic `metrics` selection for that reason, but are reported as `NA`
#' rather than dropped if you request them explicitly.
#'
#' @seealso [null_model()] for the ensembles, [summarize_graph_metrics()] for the
#'   metrics themselves.
#'
#' @examples
#' g <- igraph::sample_pa(60, power = 1.5, directed = FALSE)
#'
#' # `n_null` is small here to keep the example fast. Note the floor this puts
#' # on the attainable p-value: 1 / (10 + 1).
#' res <- metric_significance(g, metrics = c("Clustering_coefficient",
#'                                           "Modularity"),
#'                            n_null = 10, seed = 1, plot = FALSE)
#' res$result
#'
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_wrap labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom stats sd quantile
#'
#' @export
metric_significance <- function(graph,
                                null = NULL,
                                metrics = NULL,
                                n_null = 100,
                                model = c("rewire", "configuration",
                                          "erdos_renyi"),
                                weights = NULL,
                                weight_type = c("strength", "distance"),
                                seed = NULL,
                                plot = TRUE,
                                label.size = 12) {

  model <- match.arg(model)
  weight_type <- match.arg(weight_type)

  graph <- as_netkit_graph(graph)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(null)) {
    null <- null_model(graph, model = model, n = n_null,
                       shuffle_weights = !is.null(weights))
  }
  if (!is.list(null) || length(null) == 0) {
    stop("'null' must be a non-empty list of igraph objects, as returned by ",
         "null_model().", call. = FALSE)
  }

  observed <- summarize_graph_metrics(graph, weights = weights,
                                      weight_type = weight_type)

  # Null graphs carry the permuted weights under the name `weight`, whatever the
  # observed attribute was called, so that is what is requested from them.
  null_weights <- if (is.null(weights)) NULL else "weight"
  null_rows <- lapply(null, function(g) {
    suppressWarnings(
      summarize_graph_metrics(g, weights = null_weights,
                              weight_type = weight_type)
    )
  })
  null_mat <- do.call(rbind, null_rows)

  numeric_cols <- names(observed)[vapply(observed, is.numeric, logical(1))]

  if (is.null(metrics)) {
    # Drop the metrics the null model holds fixed by construction: their null
    # distribution has zero variance, so a z-score is undefined and reporting a
    # p-value for them would be meaningless rather than merely uninformative.
    varies <- vapply(numeric_cols, function(m) {
      v <- null_mat[[m]]
      sum(is.finite(v)) > 1 && stats::sd(v, na.rm = TRUE) > 0
    }, logical(1))
    metrics <- numeric_cols[varies]
    if (length(metrics) == 0) {
      stop("No metric varies across the null ensemble; nothing to test.",
           call. = FALSE)
    }
  } else {
    unknown <- setdiff(metrics, numeric_cols)
    if (length(unknown) > 0) {
      stop(sprintf(
        "Unknown or non-numeric metric(s): %s. Available: %s.",
        paste(unknown, collapse = ", "), paste(numeric_cols, collapse = ", ")
      ), call. = FALSE)
    }
  }

  rows <- lapply(metrics, function(m) {
    obs <- observed[[m]]
    nulls <- null_mat[[m]]
    nulls <- nulls[is.finite(nulls)]

    if (length(nulls) == 0 || !is.finite(obs)) {
      return(tibble::tibble(
        metric = m, observed = obs, null_mean = NA_real_, null_sd = NA_real_,
        z = NA_real_, p_empirical = NA_real_, ci_lower = NA_real_,
        ci_upper = NA_real_
      ))
    }

    mu <- mean(nulls)
    sdev <- stats::sd(nulls)

    # Two-sided, counting nulls at least as far from the null mean as the
    # observed value. The +1 / +1 convention keeps the p-value away from zero:
    # with n nulls, no result can be more significant than 1 / (n + 1).
    r <- sum(abs(nulls - mu) >= abs(obs - mu))
    p <- (r + 1) / (length(nulls) + 1)

    tibble::tibble(
      metric = m,
      observed = obs,
      null_mean = mu,
      null_sd = sdev,
      z = if (is.finite(sdev) && sdev > 0) (obs - mu) / sdev else NA_real_,
      p_empirical = p,
      ci_lower = unname(stats::quantile(nulls, 0.025)),
      ci_upper = unname(stats::quantile(nulls, 0.975))
    )
  })

  result <- do.call(rbind, rows)

  p <- NULL
  if (plot) {
    long <- do.call(rbind, lapply(metrics, function(m) {
      v <- null_mat[[m]]
      v <- v[is.finite(v)]
      if (length(v) == 0) return(NULL)
      data.frame(metric = m, value = v, stringsAsFactors = FALSE)
    }))

    if (!is.null(long) && nrow(long) > 0) {
      obs_df <- data.frame(
        metric = result$metric,
        observed = result$observed,
        stringsAsFactors = FALSE
      )
      obs_df <- obs_df[is.finite(obs_df$observed), , drop = FALSE]

      p <- ggplot2::ggplot(long, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(bins = 30, fill = "grey75", color = "white") +
        ggplot2::geom_vline(data = obs_df,
                            ggplot2::aes(xintercept = observed),
                            color = "#e41a1c", linewidth = 0.8) +
        ggplot2::facet_wrap(~ metric, scales = "free") +
        ggplot2::labs(
          x = NULL, y = "Null graphs",
          title = "Observed metrics against the null ensemble",
          subtitle = paste0(length(null), " null graphs, model '",
                            null_model_label(null, model), "'")
        ) +
        ggplot2::theme_minimal(base_size = label.size)
    }
  }

  list(
    plot = p,
    result = result,
    graph = graph,
    method = paste0(
      "Metric significance against ", length(null), " null graphs (model '",
      null_model_label(null, model), "'); ",
      describe_weights(as_netkit_weights(graph, weights, weight_type,
                                        warn_unused = FALSE)),
      "; two-sided empirical p-values with the (r+1)/(n+1) convention, so the ",
      "smallest attainable p-value is ", signif(1 / (length(null) + 1), 3)
    )
  )
}

#' Small-World Coefficients
#'
#' Computes the small-world coefficient sigma, which compares a graph's
#' clustering and path length against a degree-matched random ensemble.
#'
#' A network is "small-world" when it is much more clustered than a random graph
#' with the same degree sequence while having a comparable average path length.
#' `sigma` expresses that as a single ratio: `(C/C_rand) / (L/L_rand)`, where
#' values appreciably above 1 indicate small-world organisation.
#'
#' @param graph An `igraph` object or a data frame edge list.
#' @param n_null Integer. Number of null graphs. Default is `100`.
#' @param model Passed to [null_model()]. Default is `"rewire"`.
#' @param weights,weight_type Passed to [summarize_graph_metrics()]. See
#'   [netkit-weights].
#' @param seed Integer or `NULL`. Seed for reproducibility.
#'
#' @return A list with:
#' \describe{
#'   \item{`result`}{A one-row tibble: `sigma`, `C`, `C_rand`, `L`, `L_rand` and
#'     `n_null`.}
#'   \item{`graph`}{The input graph, unchanged.}
#'   \item{`method`}{A human-readable description.}
#' }
#'
#' @details
#' Only `sigma` is reported. The companion coefficient `omega` of Telesford et
#' al. (2011) additionally requires a *lattice* reference, and \pkg{igraph}
#' provides no degree-preserving latticisation; implementing one approximately
#' would make `omega` quietly dependent on how well that approximation worked,
#' so it is omitted rather than shipped unreliable.
#'
#' `sigma` is known to grow with network size, so it is not comparable across
#' graphs of different size -- use it to ask whether one graph is small-world,
#' not which of two is more so.
#'
#' @references
#' Humphries, M. D., & Gurney, K. (2008). Network "small-world-ness": a
#' quantitative method for determining canonical network equivalence.
#' *PLoS ONE*, 3(4), e0002051. \doi{10.1371/journal.pone.0002051}
#'
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of "small-world"
#' networks. *Nature*, 393(6684), 440-442. \doi{10.1038/30918}
#'
#' @seealso [null_model()], [metric_significance()]
#'
#' @examples
#' g <- igraph::sample_smallworld(1, 60, 4, 0.05)
#'
#' # `n_null` is small here to keep the example fast.
#' sw <- small_worldness(g, n_null = 10, seed = 1)
#' sw$result
#'
#' @importFrom tibble tibble
#'
#' @export
small_worldness <- function(graph,
                            n_null = 100,
                            model = c("rewire", "configuration",
                                      "erdos_renyi"),
                            weights = NULL,
                            weight_type = c("strength", "distance"),
                            seed = NULL) {

  model <- match.arg(model)
  weight_type <- match.arg(weight_type)

  graph <- as_netkit_graph(graph)

  if (!is.null(seed)) set.seed(seed)

  obs <- summarize_graph_metrics(graph, weights = weights,
                                 weight_type = weight_type)
  nulls <- null_model(graph, model = model, n = n_null,
                      shuffle_weights = !is.null(weights))

  null_weights <- if (is.null(weights)) NULL else "weight"
  null_stats <- do.call(rbind, lapply(nulls, function(g) {
    suppressWarnings(summarize_graph_metrics(g, weights = null_weights,
                                             weight_type = weight_type))
  }))

  C <- obs$Clustering_coefficient
  L <- obs$Average_path_length
  C_rand <- mean(null_stats$Clustering_coefficient, na.rm = TRUE)
  L_rand <- mean(null_stats$Average_path_length, na.rm = TRUE)

  # Guard every denominator rather than letting a degenerate graph produce a
  # plausible-looking sigma out of a division by zero.
  sigma <- if (!is.finite(C) || !is.finite(L) || !is.finite(C_rand) ||
               !is.finite(L_rand) || C_rand == 0 || L_rand == 0 || L == 0) {
    NA_real_
  } else {
    (C / C_rand) / (L / L_rand)
  }

  list(
    result = tibble::tibble(
      sigma = sigma,
      C = C, C_rand = C_rand,
      L = L, L_rand = L_rand,
      n_null = length(nulls)
    ),
    graph = graph,
    method = paste0(
      "Small-world sigma = (C/C_rand)/(L/L_rand) against ", length(nulls),
      " null graphs (model '", model, "'); ",
      describe_weights(as_netkit_weights(graph, weights, weight_type,
                                         warn_unused = FALSE))
    )
  )
}

#' Report which null model an ensemble came from
#'
#' Internal. [null_model()] records its model in an attribute, but a caller may
#' pass any list of graphs, so fall back to the `model` argument rather than
#' reporting NULL.
#'
#' @param null The ensemble.
#' @param model The `model` argument of the calling function.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
null_model_label <- function(null, model) {
  lab <- attr(null, "model")
  if (is.null(lab)) model else lab
}
