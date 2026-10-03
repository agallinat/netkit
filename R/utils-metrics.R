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
