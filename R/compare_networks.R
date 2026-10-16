#' Compare Two Networks
#'
#' This function compares two networks using summary metrics, degree distributions,
#' and topological similarity measures. It overlays the complementary cumulative
#' frequency distributions (CCDFs) of degree and returns a combined report.
#'
#' @param graph1 An igraph object or data.frame (edge list).
#' @param graph2 An igraph object or data.frame (edge list).
#' @param remove_singles Logical; remove single nodes before analysis.
#' @param label.size Labels' size in the CCDF plot.
#' @param show_PL Logical; whether to fit and display power law exponents.
#' @param PL_exponents Vector; power-law slopes to show.
#' @param colors Optional vector of colors for the CCDF plot.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{plot}}{A \code{ggplot2} object overlaying both CCDF curves.}
#'   \item{\code{global_topology}}{A data frame with one row per input graph, as
#'     produced by [summarize_graph_metrics()].}
#'   \item{\code{similarity}}{A one-row data frame with three distinct set
#'     statistics: \code{jaccard_similarity}, the Jaccard index of the edge sets
#'     (shared edges over their union); \code{node_overlap}, the Jaccard index of
#'     the vertex sets; and \code{edge_overlap}, the overlap coefficient of the
#'     edge sets (shared edges over the \emph{smaller} of the two edge sets).
#'     The last is the informative companion to Jaccard when the two networks
#'     differ greatly in size -- a small network nested inside a large one scores
#'     near 1 on overlap and near 0 on Jaccard. Each is \code{NaN} where its
#'     denominator is empty.}
#'   \item{\code{edge_enrichment}}{A one-row data frame placing
#'     \code{shared_edges} against the number expected if the two edge sets were
#'     drawn independently. See the "Shared edges against chance" section.}
#'   \item{\code{ks_test}}{The Kolmogorov-Smirnov test comparing the two degree
#'     distributions, as returned by [stats::ks.test()].}
#' }
#'
#' @section Shared edges against chance:
#'
#' A count of shared edges is not interpretable on its own: two large dense
#' networks share many edges by arithmetic alone, and the question is always
#' whether they share more than that. `edge_enrichment` answers it by treating
#' each network's edge set as a draw from the pool of vertex pairs the two
#' networks could possibly agree on:
#'
#' \describe{
#'   \item{\code{n_shared_nodes}}{Vertices present in both networks. An edge can
#'     only be shared if both its endpoints are, so this is the vertex set the
#'     test is conducted on -- not the union, which would include pairs neither
#'     network could ever have reported.}
#'   \item{\code{n_possible_edges}}{\eqn{\binom{n}{2}} over those vertices: the
#'     size of the pool.}
#'   \item{\code{edges1}, \code{edges2}}{Distinct non-loop edges of each network
#'     lying inside that vertex set -- the two draws. These are generally smaller
#'     than the `Edges` column of `global_topology`, which counts every edge of
#'     each network whether or not its endpoints are shared.}
#'   \item{\code{shared_edges}}{Edges drawn by both.}
#'   \item{\code{expected_shared_edges}}{\eqn{m_1 m_2 / \binom{n}{2}}, the
#'     expectation under independence.}
#'   \item{\code{fold_enrichment}}{Observed over expected. Above 1 means the two
#'     networks agree more than chance, below 1 less.}
#'   \item{\code{p_value}}{A one-sided hypergeometric (Fisher) upper-tail
#'     probability of seeing at least this many shared edges, from
#'     [stats::phyper()]. One-sided: it tests enrichment, not depletion.}
#' }
#'
#' Direction and multiplicity are collapsed before counting, matching the
#' treatment in `similarity` -- an edge is an unordered pair of vertex
#' identifiers, counted once. Self-loops are excluded, since
#' \eqn{\binom{n}{2}} does not contain them.
#'
#' The independence assumption is the usual one for this test and is only
#' approximately true of networks: both edge sets are constrained by their degree
#' sequences, and two networks over the same vertices tend to share hubs. Read
#' `fold_enrichment` as a scale for the overlap rather than as a calibrated
#' significance test; for a calibrated one, permute within the degree sequence
#' using [null_model()].
#'
#' Where the pool is empty -- fewer than two shared vertices -- the expectation
#' is `NaN` and `p_value` is `NA`, rather than a number computed from a
#' denominator of zero.
#'
#' @examples
#' g1 <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
#' g2 <- igraph::sample_gnp(80, 0.05, directed = FALSE)
#' res <- compare_networks(g1, g2)
#' res$global_topology
#' res$ks_test
#'
#' # Two unrelated networks over the same vertices: the shared edges are about
#' # as many as chance gives, so fold_enrichment sits near 1.
#' res$edge_enrichment
#'
#' @importFrom dplyr bind_rows
#' @importFrom igraph is_igraph degree induced_subgraph graph_from_data_frame V E
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual labs coord_cartesian theme_minimal scale_y_log10 scale_x_continuous
#' @importFrom rlang sym
#' @importFrom scales trans_breaks trans_format math_format label_math hue_pal
#' @importFrom stats setNames ks.test
#' @param weights1,weights2 Optional edge weights for each graph: `NULL`
#'   (default) to ignore them, the name of an edge attribute, or a numeric
#'   vector. Supplied separately because the two graphs need not carry the same
#'   attribute. See [netkit-weights].
#' @param weight_type Either `"strength"` (default) or `"distance"`, applied to
#'   both graphs. See [netkit-weights].
#'
#' @inheritSection netkit-weights Edge weights
#'
#' @importFrom graphics par
#'
#' @export
compare_networks <- function(graph1, graph2,
                             remove_singles = FALSE,
                             show_PL = TRUE,
                             PL_exponents = c(2, 3),
                             colors = c("#e41a1c", "#000831", "#9c52f2", "#b8b8ff"),
                             label.size = 12,
                             weights1 = NULL,
                             weights2 = NULL,
                             weight_type = c("strength", "distance")
) {

  weight_type <- match.arg(weight_type)

  # --- Validate input ---
  graph1 <- as_netkit_graph(graph1, arg = "graph1")
  graph2 <- as_netkit_graph(graph2, arg = "graph2")

  # --- Handle single nodes ---
  if (remove_singles) {
    degree_g1 <- igraph::degree(graph1)
    graph1 <- induced_subgraph(graph1, vids = which(degree_g1 > 0))
    degree_g2 <- igraph::degree(graph2)
    graph2 <- induced_subgraph(graph2, vids = which(degree_g2 > 0))
    if (0 %in% c(degree_g1, degree_g2)) {
      message("Single nodes excluded from the analysis. ",
              "Set 'remove_singles' to FALSE to include all nodes.")
    }
  }

  # --- Compute basic metrics ---
  metrics1 <- summarize_graph_metrics(graph1, weights = weights1,
                                      weight_type = weight_type)
  metrics2 <- summarize_graph_metrics(graph2, weights = weights2,
                                      weight_type = weight_type)

  # --- Degree distributions ---
  deg1 <- igraph::degree(graph1)
  deg2 <- igraph::degree(graph2)

  # --- KS Test ---
  ks <- suppressWarnings(stats::ks.test(deg1, deg2))

  # --- Similarity ---
  nodes1 <- igraph::vertex_attr(graph1, "name")
  nodes2 <- igraph::vertex_attr(graph2, "name")
  shared_nodes <- intersect(nodes1, nodes2)
  node_overlap <- length(shared_nodes) / length(union(nodes1, nodes2))

  # Vertex identifiers and canonical edge keys come from the same place, so the
  # edge sets and the vertex universe used for the enrichment test below cannot
  # drift apart. `names = FALSE` then an explicit lookup, rather than letting
  # as_edgelist() choose: for an unnamed graph it returns doubles, and
  # paste()-ing those goes through scientific notation above 1e5, which would
  # key the same vertex differently from as.character(seq_len(.)).
  ids1 <- graph_vertex_ids(graph1)
  ids2 <- graph_vertex_ids(graph2)
  edge_set1 <- edge_keys(ids1, igraph::as_edgelist(graph1, names = FALSE))
  edge_set2 <- edge_keys(ids2, igraph::as_edgelist(graph2, names = FALSE))
  shared_edges <- length(intersect(edge_set1, edge_set2))

  # Jaccard and the overlap coefficient answer different questions, and
  # `edge_overlap` used to be a duplicate of `jaccard_similarity` -- the identical
  # expression under a second name. Dividing by the smaller edge set instead makes
  # it the overlap (Szymkiewicz-Simpson) coefficient, which is the informative
  # companion when the two networks differ greatly in size: a small curated
  # network nested inside a large screen scores near 1 here and near 0 on Jaccard.
  union_edges <- length(union(edge_set1, edge_set2))
  min_edges <- min(length(unique(edge_set1)), length(unique(edge_set2)))

  jaccard_sim <- if (union_edges == 0) NaN else shared_edges / union_edges
  edge_overlap <- if (min_edges == 0) NaN else shared_edges / min_edges

  similarity <- data.frame(jaccard_similarity = jaccard_sim,
                           node_overlap = node_overlap,
                           edge_overlap = edge_overlap)

  # --- Shared edges against chance ---
  #
  # The three statistics above are all ratios of observed quantities, and none of
  # them says whether the overlap is larger than arithmetic alone would give:
  # two dense networks over the same vertices share edges whether or not they
  # agree about anything. Reported as its own element rather than as further
  # columns of `similarity`, which is documented as a frame of proportions and
  # is range-checked as one.
  edge_enrichment <- shared_edge_enrichment(ids1, ids2, edge_set1, edge_set2)

  # Compute ccdfs. rep() rather than a scalar: compute_ccdf() returns zero rows
  # for an edgeless graph, and recycling a length-1 value into a zero-row frame
  # is an error ("replacement has 1 row, data has 0").
  ccdf1 <- compute_ccdf(graph1)
  ccdf1$graph <- rep("Graph 1", nrow(ccdf1))

  ccdf2 <- compute_ccdf(graph2)
  ccdf2$graph <- rep("Graph 2", nrow(ccdf2))

  ccdf_combined <- rbind(ccdf1, ccdf2)

  # Ensure 'colors' has enough entries
  n_lines <- if (show_PL) length(PL_exponents) + 2 else 2

  if (length(colors) < n_lines) {
    colors <- hue_pal(l = 65, c = 100)(n_lines)
  }

  colors_vec <- c("Graph 1" = colors[1], "Graph 2" = colors[2])

  if (show_PL) {
    for (i in seq_along(PL_exponents)) {
      gamma <- PL_exponents[i]
      ccdf_combined[[paste0("PL", gamma)]] <- ccdf_combined$degree^(-gamma)
      colors_vec <- c(colors_vec, setNames(colors[i + 2], paste0("gamma = ", gamma)))
    }
  }

  # Plot overlay
  p_combined <- ggplot(ccdf_combined, aes(x = degree, y = ccdf, color = graph)) +
    geom_line(aes(y=ccdf), linewidth = 0.7)

  if (show_PL) {
    for (gamma in PL_exponents) {
      col_name <- paste0("PL", gamma)
      label <- paste0("gamma = ", gamma)
      # `!!` forces col_name and label at aes() construction time; see plot_CCDF().
      p_combined <- p_combined + geom_line(aes(y = !!sym(col_name), color = !!label),
                                           linetype = "dashed", linewidth = 0.5)
    }
  }

  p_combined <- p_combined +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^floor(x)),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      transform = "log2") +
    labs(x = "Degree, k", y = "Pr(K > k)", color = "Graph") +
    scale_color_manual(values = colors_vec) +
    theme_minimal(base_size = label.size)

  return(as_netkit_result(list(
    plot = p_combined,
    global_topology = rbind(metrics1, metrics2),
    similarity = similarity,
    edge_enrichment = edge_enrichment,
    ks_test = ks
  ), "compare_networks"))
}

#' Vertex identifiers as the edge list reports them
#'
#' Internal. Returns `V(graph)$name` when the graph is named, and positional
#' indices as character otherwise, which is what `igraph::as_edgelist()` falls
#' back to. Used so that edge keys and the vertex universe they are tested
#' against are built from one source.
#'
#' @param graph An `igraph` object.
#'
#' @return A character vector of length `igraph::vcount(graph)`.
#'
#' @keywords internal
#' @noRd
graph_vertex_ids <- function(graph) {
  nm <- igraph::vertex_attr(graph, "name")
  if (is.null(nm)) as.character(seq_len(igraph::vcount(graph))) else as.character(nm)
}

#' Canonical undirected edge keys
#'
#' Internal. Collapses each edge to a `"lo|hi"` string over sorted endpoint
#' identifiers, so that an edge is an unordered pair and direction is ignored.
#' Vectorized through `pmin()`/`pmax()` rather than `apply(., 1, sort)`, which
#' called an R closure per edge.
#'
#' @param ids Character vertex identifiers, as from `graph_vertex_ids()`.
#' @param el A two-column edge list of vertex *indices*.
#'
#' @return A character vector with one key per edge, with duplicates retained.
#'
#' @keywords internal
#' @noRd
edge_keys <- function(ids, el) {
  if (nrow(el) == 0) return(character(0))
  a <- ids[el[, 1]]
  b <- ids[el[, 2]]
  paste(pmin(a, b), pmax(a, b), sep = "|")
}

#' Shared edges against their random expectation
#'
#' Internal. Places the observed edge overlap against the count expected if the
#' two edge sets were independent draws from the pairs the two networks could
#' both have reported.
#'
#' The universe is \eqn{\binom{n}{2}} over the *shared* vertices, and the two
#' draw sizes are each network's distinct non-loop edges **within that vertex
#' set**. Deriving the counts and the denominator from the same restriction is
#' the point: taking \eqn{m_1} from `ecount()` while taking \eqn{n} from the
#' intersection mixes two different universes and inflates the expectation by
#' however many edges reach a vertex the other network never had -- the same
#' count-and-denominator mismatch that produced the `calculate_roles()` and
#' `robustness_analysis()` defects.
#'
#' @param ids1,ids2 Vertex identifiers of the two graphs.
#' @param keys1,keys2 Canonical edge keys, as from `edge_keys()`.
#'
#' @return A one-row data frame; see the "Shared edges against chance" section
#'   of [compare_networks()].
#'
#' @keywords internal
#' @noRd
shared_edge_enrichment <- function(ids1, ids2, keys1, keys2) {

  shared_ids <- intersect(ids1, ids2)
  n_shared <- length(shared_ids)
  n_possible <- if (n_shared < 2) 0 else choose(n_shared, 2)

  # An edge qualifies only if both endpoints are shared and it is not a loop.
  # strsplit() on the key rather than carrying the endpoints along: the key is
  # the canonical form, and re-deriving the endpoints from the edge list would
  # risk the two disagreeing about multi-edges.
  restrict <- function(keys) {
    if (length(keys) == 0) return(character(0))
    ends <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
    keep <- ends[, 1] != ends[, 2] &
      ends[, 1] %in% shared_ids & ends[, 2] %in% shared_ids
    unique(keys[keep])
  }

  e1 <- restrict(keys1)
  e2 <- restrict(keys2)
  m1 <- length(e1)
  m2 <- length(e2)
  shared <- length(intersect(e1, e2))

  expected <- if (n_possible == 0) NaN else m1 * m2 / n_possible
  fold <- if (is.nan(expected) || expected == 0) NaN else shared / expected

  # Upper tail: P(X >= shared) = 1 - P(X <= shared - 1). One-sided on purpose --
  # the question asked of two networks is whether they agree more than chance.
  p <- if (n_possible == 0) {
    NA_real_
  } else {
    stats::phyper(shared - 1, m1, n_possible - m1, m2, lower.tail = FALSE)
  }

  data.frame(
    n_shared_nodes = n_shared,
    n_possible_edges = n_possible,
    edges1 = m1,
    edges2 = m2,
    shared_edges = shared,
    expected_shared_edges = expected,
    fold_enrichment = fold,
    p_value = p
  )
}
