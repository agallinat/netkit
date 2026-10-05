#' Printing netkit objects
#'
#' netkit's analysis functions return a list holding a `result` table, the
#' annotated `graph`, a `method` string describing what was computed and, where
#' there is one, a `plot`. Printing such a list with R's default method dumps
#' the whole graph and every row of the table, which for a 100-graph
#' [null_model()] ensemble runs to well over a thousand lines. These methods
#' print a summary instead.
#'
#' Only the printing changes. The objects are still plain lists: `x$result`,
#' `x$graph`, `x$plot` and `x[[i]]` behave exactly as they did, and
#' `unclass(x)` recovers the default printing.
#'
#' @param x A netkit object: an analysis result, a [null_model()] ensemble or a
#'   [prepare_diffusion()] kernel.
#' @param n Number of rows of `x$result` to preview. Default 5.
#' @param ... Ignored, present for consistency with the generic.
#'
#' @return `x`, invisibly. Called for the side effect of printing.
#'
#' @examples
#' g <- igraph::sample_gnp(40, 0.1, directed = FALSE)
#' igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
#'
#' find_hubs(g, plot = FALSE)
#'
#' # The elements are unchanged -- only the printing is summarized.
#' hubs <- find_hubs(g, plot = FALSE)
#' head(hubs$result, 3)
#'
#' # n = 20 rather than the default 100 null graphs, to keep the example quick.
#' null_model(g, n = 20, seed = 1)
#'
#' prepare_diffusion(g, method = "rwr")
#'
#' @name netkit-print
NULL

#' Attach the netkit result class
#'
#' Internal. Classing the returned list is what lets [print.netkit_result()]
#' summarize it. The class vector carries a per-function subclass
#' (`"netkit_find_hubs"` and so on) so that a caller can dispatch on a specific
#' result, and `"list"` last so that anything testing for a plain list still
#' sees one.
#'
#' @param x The list to class.
#' @param fn Name of the function producing it, used in the printed header.
#'
#' @return `x`, classed.
#'
#' @keywords internal
#' @noRd
as_netkit_result <- function(x, fn) {
  attr(x, "netkit_fn") <- fn
  class(x) <- c(paste0("netkit_", fn), "netkit_result", "list")
  x
}

#' Shorten a string for one-line display
#'
#' Internal.
#'
#' @param x A single string.
#' @param max Maximum number of characters to keep.
#'
#' @return A single string, suffixed with an ellipsis if it was shortened.
#'
#' @keywords internal
#' @noRd
truncate_chr <- function(x, max = 48) {
  if (nchar(x) <= max) return(x)
  paste0(substr(x, 1, max - 3), "...")
}

#' One-line description of a list element
#'
#' Internal. Used by [print.netkit_result()] for every element it does not
#' print in full, so that a user can see what is there without the console
#' filling with an adjacency matrix.
#'
#' @param el The element to describe.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
describe_slot <- function(el) {
  if (is.null(el)) {
    return("NULL")
  }
  # ggplot2 >= 4.0 reports its S7 class as "ggplot2::ggplot"; the namespace adds
  # nothing here.
  cls <- sub(".*::", "", class(el)[1])
  if (igraph::is_igraph(el)) {
    return(sprintf(
      "<igraph> %d nodes, %d edges, %s | vertex attrs: %s",
      igraph::vcount(el), igraph::ecount(el),
      if (igraph::is_directed(el)) "directed" else "undirected",
      truncate_chr(paste(igraph::vertex_attr_names(el), collapse = ", "))
    ))
  }
  if (inherits(el, "data.frame")) {
    return(sprintf("<%s> %d x %d | %s", cls, nrow(el), ncol(el),
                   truncate_chr(paste(names(el), collapse = ", "))))
  }
  # ggExtra::ggMarginal() returns an assembled gtable rather than a ggplot, so
  # both spellings have to be recognized here.
  if (inherits(el, c("ggplot", "ggExtraPlot", "gtable", "grob"))) {
    return(sprintf("<%s> print(x$plot) to draw", cls))
  }
  if (is.atomic(el)) {
    if (length(el) == 1) {
      return(truncate_chr(format(el)))
    }
    return(sprintf("<%s [%d]> %s", typeof(el), length(el),
                   truncate_chr(paste(format(utils::head(el, 5)),
                                      collapse = ", "))))
  }
  if (is.list(el)) {
    nms <- names(el)
    return(sprintf("<%s [%d]>%s", cls, length(el),
                   if (is.null(nms)) "" else
                     paste0(" ", truncate_chr(paste(nms, collapse = ", ")))))
  }
  paste0("<", cls, ">")
}

#' @rdname netkit-print
#' @export
print.netkit_result <- function(x, n = 5, ...) {
  fn <- attr(x, "netkit_fn")
  cat("<netkit result: ", if (is.null(fn)) "unknown" else paste0(fn, "()"),
      ">\n", sep = "")

  # `method` records the thresholds actually used, which is the thing most
  # worth seeing without asking for it, so it goes in the header rather than
  # the element list.
  shown <- character(0)
  if (!is.null(x$method) && is.character(x$method)) {
    cat(paste(strwrap(paste(x$method, collapse = " "),
                      width = max(40, getOption("width", 80) - 2),
                      initial = "Method: ", prefix = "  "),
              collapse = "\n"), "\n", sep = "")
    shown <- "method"
  }

  if (inherits(x$result, "data.frame")) {
    cat("\n$result\n")
    print(tibble::as_tibble(x$result), n = n)
    shown <- c(shown, "result")
  }

  rest <- setdiff(names(x), shown)
  if (length(rest) > 0) {
    if (length(shown) > 0) cat("\n")
    pad <- max(nchar(rest))
    for (nm in rest) {
      cat(sprintf("$%-*s  %s\n", pad, nm, describe_slot(x[[nm]])))
    }
  }

  invisible(x)
}

#' @rdname netkit-print
#' @export
print.netkit_null <- function(x, ...) {
  model <- attr(x, "model")
  cat(sprintf("<netkit null ensemble: %d graph%s, model '%s'>\n",
              length(x), if (length(x) == 1) "" else "s",
              if (is.null(model)) "unknown" else model))

  if (length(x) > 0 && igraph::is_igraph(x[[1]])) {
    cat(sprintf("Each graph has %d nodes and %d edges.\n",
                igraph::vcount(x[[1]]), igraph::ecount(x[[1]])))
  }
  cat("Use x[[i]] for one graph, or pass the ensemble to",
      "metric_significance(null = x).\n")

  invisible(x)
}

#' @param i Indices of the null graphs to keep.
#'
#' @rdname netkit-print
#' @export
`[.netkit_null` <- function(x, i) {
  # Without this, `nulls[1:10]` drops the class and the model attribute, and
  # printing the subset dumps every graph again -- the thing print.netkit_null()
  # exists to prevent. A subset of an ensemble is still an ensemble.
  structure(unclass(x)[i], class = c("netkit_null", "list"),
            model = attr(x, "model"))
}

#' @rdname netkit-print
#' @export
print.netkit_kernel <- function(x, ...) {
  cat(sprintf("<netkit diffusion kernel: method '%s'>\n", x$method))
  if (!is.null(x$L)) {
    cat(sprintf("%d x %d Laplacian, %s edges.\n", nrow(x$L), ncol(x$L),
                if (identical(x$weights_key, "unweighted")) "unweighted" else
                  "weighted"))
  }
  cat("Pass to network_diffusion(precompute = x) to skip rebuilding it.\n")

  invisible(x)
}
