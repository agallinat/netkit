#' Participation-coefficient boundaries for the Guimera-Amaral roles
#'
#' Internal. The seven connectivity roles are delimited by a within-module
#' z-score threshold (the `hub_z` argument of [calculate_roles()]) and by four
#' participation-coefficient boundaries, published in Guimera & Amaral (2005).
#'
#' These numbers are defined here once because they previously appeared in three
#' places inside `calculate_roles()` -- the classifier, the `roles_definitions`
#' table handed back to the caller, and the shaded bands of the diagnostic plot --
#' and had drifted apart. The R5/R6 boundary was `0.30` in the classifier but
#' `0.25` in both the table and the plot, so the shaded "provincial hub" region
#' did not match the classification it was illustrating, and the table documented
#' a threshold the code did not use.
#'
#' @return A named numeric vector of length four, in increasing order within each
#'   of the two z-score strata:
#'   \describe{
#'     \item{`R1_R2`}{Non-hub: ultra-peripheral below, peripheral above (0.05).}
#'     \item{`R2_R3`}{Non-hub: peripheral below, connector above (0.62).}
#'     \item{`R3_R4`}{Non-hub: connector below, kinless above (0.80).}
#'     \item{`R5_R6`}{Hub: provincial below, connector above (0.30).}
#'     \item{`R6_R7`}{Hub: connector below, kinless above (0.75).}
#'   }
#'
#' @keywords internal
#' @noRd
netkit_role_p <- function() {
  c(R1_R2 = 0.05, R2_R3 = 0.62, R3_R4 = 0.80, R5_R6 = 0.30, R6_R7 = 0.75)
}

#' Validate and merge user-supplied role thresholds
#'
#' Internal. Resolves the `thresholds` argument of [calculate_roles()] against
#' the published defaults, allowing a partial override while keeping the result a
#' complete, ordered, in-range set.
#'
#' @param thresholds `NULL` for the published defaults, or a named numeric vector
#'   overriding one or more of them. Names must be a subset of those returned by
#'   `netkit_role_p()`.
#'
#' @return The complete named numeric vector of boundaries.
#'
#' @keywords internal
#' @noRd
as_role_thresholds <- function(thresholds = NULL) {

  defaults <- netkit_role_p()

  if (is.null(thresholds)) {
    return(defaults)
  }

  if (!is.numeric(thresholds) || is.null(names(thresholds))) {
    stop("'thresholds' must be a named numeric vector, or NULL for the ",
         "Guimera-Amaral (2005) defaults.", call. = FALSE)
  }

  unknown <- setdiff(names(thresholds), names(defaults))
  if (length(unknown) > 0) {
    stop(sprintf(
      "Unknown 'thresholds' name(s): %s. Valid names are: %s.",
      paste(unknown, collapse = ", "), paste(names(defaults), collapse = ", ")
    ), call. = FALSE)
  }

  if (anyNA(thresholds) || any(thresholds < 0) || any(thresholds > 1)) {
    stop("'thresholds' values must be non-missing and within [0, 1]: the ",
         "participation coefficient is bounded by 0 and 1.", call. = FALSE)
  }

  out <- defaults
  out[names(thresholds)] <- as.numeric(thresholds)

  # The two z-score strata are classified by independent if/else ladders, so each
  # must be monotonically increasing on its own for every role to be reachable.
  if (!all(diff(out[c("R1_R2", "R2_R3", "R3_R4")]) > 0)) {
    stop("'thresholds' R1_R2 < R2_R3 < R3_R4 must hold, otherwise a non-hub ",
         "role is unreachable.", call. = FALSE)
  }
  if (!all(diff(out[c("R5_R6", "R6_R7")]) > 0)) {
    stop("'thresholds' R5_R6 < R6_R7 must hold, otherwise a hub role is ",
         "unreachable.", call. = FALSE)
  }

  out
}
