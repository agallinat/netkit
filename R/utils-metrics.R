#' Scale a metric series to its own maximum
#'
#' Internal helper. Robustness curves are plotted and integrated on a 0-1 scale by
#' dividing each metric by its maximum over the removal sequence. Dividing directly
#' is unsafe: if the maximum is zero (for example global efficiency on a graph with
#' no edges, which is zero at every step) every value becomes `NaN`, which then
#' propagates silently into [pracma::trapz()] and produces a `NaN` AUC.
#'
#' @param x Numeric vector of metric values, possibly containing `NA`/`NaN`.
#'
#' @return `x` divided by its maximum, or `x` unchanged when that maximum is zero or
#'   not finite.
#'
#' @keywords internal
#' @noRd
normalize_metric <- function(x) {
  m <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(m) || m == 0) {
    return(x)
  }
  x / m
}

#' Evaluate an expression under a fixed seed, leaving the caller's stream alone
#'
#' Internal helper. Some metrics are computed by a randomized algorithm --
#' Louvain visits vertices in a random order -- which makes them a random
#' variable rather than a property of the graph: repeated
#' `summarize_graph_metrics()` calls on one graph returned modularities spanning
#' 0.356 to 0.366, so a reported "the modularity is 0.36" was not reproducible
#' and two summaries could not be compared.
#'
#' Seeding globally would fix that at the cost of silently resetting the
#' caller's random stream, which is the defect `robustness_analysis()` and
#' `null_model()` were both corrected for. This saves `.Random.seed`, seeds,
#' evaluates, and restores -- so the metric is deterministic and nothing outside
#' can tell the function used the RNG at all.
#'
#' @param seed Integer seed, or `NULL` to evaluate under the ambient stream
#'   without touching it.
#' @param expr Expression to evaluate. Lazily passed and forced after seeding.
#'
#' @return The value of `expr`.
#'
#' @keywords internal
#' @noRd
with_local_seed <- function(seed, expr) {

  if (is.null(seed)) {
    return(force(expr))
  }

  # R only materializes .Random.seed once the RNG has been used; if it is absent
  # the honest restoration is to remove it again, not to leave ours behind.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }

  set.seed(seed)
  force(expr)
}
