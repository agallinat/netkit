#' Standardize a node metric, reporting a degenerate distribution
#'
#' Internal. `scale()` divides by the standard deviation, which is zero whenever
#' every vertex carries the same value -- a regular graph, a ring, a complete
#' graph, or any graph whose betweenness is uniformly zero. It then returns
#' `NaN` for every vertex, and a comparison against `NaN` is `NA`, so
#' `find_hubs()` and `find_bottlenecks()` used to return an `is_hub` /
#' `is_bottleneck` column that was `NA` throughout and write that onto the graph
#' as a vertex attribute. Downstream, `table()` dropped it and a filter returned
#' nothing, so the graph looked like it had no hubs rather than like a graph the
#' question does not apply to.
#'
#' The z-scores stay `NA`, which is what an undefined standardization should
#' report -- the same convention `metric_significance()` uses for a zero-variance
#' null. The *classification* does not: "is this vertex unusually central" has a
#' defined answer when no vertex differs from any other, and the answer is no.
#' The caller substitutes `FALSE` and the flag returned here drives the warning.
#'
#' @param x A numeric vector of per-vertex metric values.
#'
#' @return A list with `z` (numeric, `NA` throughout when degenerate) and
#'   `degenerate` (logical).
#'
#' @keywords internal
#' @noRd
standardize_node_metric <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(list(z = rep(NA_real_, length(x)), degenerate = TRUE))
  }
  list(z = as.numeric(scale(x)[, 1]), degenerate = FALSE)
}

#' Warn when a two-axis threshold pair selects no vertices
#'
#' Internal. Shared by `find_hubs()` and `find_bottlenecks()`, which select on
#' the conjunction of a degree criterion and a betweenness criterion and
#' previously returned an empty selection in silence. An empty selection has
#' three quite different causes, and the remedy differs for each, so the warning
#' diagnoses which one occurred rather than merely reporting the count:
#'
#' * A metric with zero variance, where the standardization is undefined.
#' * A threshold no vertex can reach. This is what the `"zscore"` defaults do on
#'   a dense graph: the degree distribution is near-Gaussian, so the largest
#'   attainable z-score is around 2.5 for a few hundred vertices and a threshold
#'   of 3 is unreachable rather than unmet. It cannot happen under `"quantile"`,
#'   where the cutoff is by construction a value some vertex attains.
#' * Two criteria each met by some vertex, but never by the same one. That is a
#'   statement about the graph and not about the arguments, so it is reported
#'   without suggesting a change.
#'
#' @param fn The calling function's name, for the message.
#' @param label Plural noun for what was not found, e.g. `"hub nodes"`.
#' @param method The thresholding method in use.
#' @param axes A list with one element per criterion, each a list of `metric`
#'   (the argument name), `label` (human-readable), `direction` (`">="` or
#'   `"<="`), `threshold`, `attained` (the most favourable value any vertex
#'   reached: the maximum for `">="`, the minimum for `"<="`) and `degenerate`.
#'
#' @return Invisibly `NULL`; called for the warning.
#'
#' @keywords internal
#' @noRd
warn_empty_selection <- function(fn, label, method, axes) {

  fmt <- function(x) format(signif(unname(x), 3))

  degenerate <- vapply(axes, function(a) isTRUE(a$degenerate), logical(1))

  if (any(degenerate)) {
    flat <- vapply(axes[degenerate], function(a) a$label, character(1))
    warning(sprintf(
      paste0("%s() found no %s: %s %s the same value at every vertex, so the ",
             "standardized metric is undefined and no vertex can be unusual. ",
             "This is a property of the graph -- in a regular graph no vertex ",
             "is better connected than another. The standardized columns are ",
             "NA for that reason; the classification is FALSE throughout."),
      fn, label, paste(flat, collapse = " and "),
      if (length(flat) > 1) "take" else "takes"
    ), call. = FALSE)
    return(invisible(NULL))
  }

  unreachable <- vapply(axes, function(a) {
    if (a$direction == ">=") a$attained < a$threshold else a$attained > a$threshold
  }, logical(1))

  if (any(unreachable)) {
    clauses <- vapply(axes[unreachable], function(a) {
      sprintf("no vertex reaches %s = %s (the %s %s is %s)",
              a$metric, fmt(a$threshold),
              if (a$direction == ">=") "highest" else "lowest",
              a$label, fmt(a$attained))
    }, character(1))

    advice <- if (method == "zscore") {
      paste0(" A z-score threshold is only reachable if the metric has a ",
             "heavy enough tail: on a dense graph the degree distribution is ",
             "near-Gaussian and the largest attainable z-score is around 2 to ",
             "3, so a threshold above that is unreachable rather than unmet. ",
             "Use method = \"quantile\" for percentile cutoffs, or lower the ",
             "threshold.")
    } else {
      " Lower the corresponding quantile to widen the selection."
    }

    warning(sprintf("%s() found no %s: %s.%s", fn, label,
                    paste(clauses, collapse = ", and "), advice),
            call. = FALSE)
    return(invisible(NULL))
  }

  # Each criterion is met by someone, just never by the same vertex. Nothing to
  # advise: this is the honest answer to the question that was asked.
  warning(sprintf(
    paste0("%s() found no %s: some vertices pass %s = %s and some pass ",
           "%s = %s, but no vertex passes both. The two criteria are met by ",
           "disjoint sets of vertices in this graph."),
    fn, label,
    axes[[1]]$metric, fmt(axes[[1]]$threshold),
    axes[[2]]$metric, fmt(axes[[2]]$threshold)
  ), call. = FALSE)

  invisible(NULL)
}
