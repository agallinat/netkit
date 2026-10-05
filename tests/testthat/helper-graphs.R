# Shared, deterministic fixtures for the test suite.
#
# Every graph here is small and built under a fixed seed so that tests are fast
# and reproducible. netkit's community detection, random node removal and
# permutation nulls are all stochastic, so tests that touch them must either set
# a seed or assert only on structure rather than on exact values.

# A scale-free-ish graph with named vertices and a numeric vertex attribute.
# Scale-free degree distributions give find_hubs()/find_bottlenecks() something
# to actually find.
test_graph <- function(n = 60, seed = 123) {
  set.seed(seed)
  g <- igraph::sample_pa(n, power = 1.5, directed = FALSE)
  igraph::V(g)$name <- paste0("n", seq_len(igraph::vcount(g)))
  igraph::V(g)$score <- stats::rnorm(igraph::vcount(g))
  g
}

# A denser Erdos-Renyi graph, used as the second network in comparisons.
test_graph_gnp <- function(n = 50, p = 0.08, seed = 456) {
  set.seed(seed)
  g <- igraph::sample_gnp(n, p, directed = FALSE)
  igraph::V(g)$name <- paste0("m", seq_len(igraph::vcount(g)))
  g
}

# A directed graph, for the undirected-coercion paths.
test_graph_directed <- function(n = 40, p = 0.06, seed = 789) {
  set.seed(seed)
  g <- igraph::sample_gnp(n, p, directed = TRUE)
  igraph::V(g)$name <- paste0("d", seq_len(igraph::vcount(g)))
  g
}

# A 5-node ring: every global metric is known by hand, so it pins
# summarize_graph_metrics() to exact values rather than to its own output.
ring_graph <- function() {
  g <- igraph::make_ring(5)
  igraph::V(g)$name <- letters[1:5]
  g
}

# Edge-list form of a graph, for the data.frame input path.
as_edge_df <- function(g) {
  igraph::as_data_frame(g, what = "edges")
}

# A weighted graph. Weights are a *strength* (confidence-score style, higher =
# stronger), which is netkit's default interpretation, and they are deliberately
# spread over two orders of magnitude so that ignoring them is detectable.
test_graph_weighted <- function(n = 60, seed = 321) {
  set.seed(seed)
  g <- igraph::sample_gnp(n, 0.08, directed = FALSE)
  igraph::V(g)$name <- paste0("w", seq_len(igraph::vcount(g)))
  igraph::E(g)$weight <- stats::runif(igraph::ecount(g), 0.01, 1)
  g
}

# The same graph with every weight equal to 1. Weighted and unweighted results
# must agree here, which is the single property that pins the whole weight
# contract: it holds for every metric at once and needs no reference values.
test_graph_unit_weights <- function(n = 60, seed = 321) {
  g <- test_graph_weighted(n = n, seed = seed)
  igraph::E(g)$weight <- rep(1, igraph::ecount(g))
  g
}

# Two disjoint cliques joined by nothing: a deliberately disconnected graph for
# the Steiner/closeness/spinglass paths, where "no path exists" is the case that
# matters and the component structure is known by hand.
test_graph_disconnected <- function() {
  g <- igraph::disjoint_union(igraph::make_full_graph(5), igraph::make_full_graph(4))
  igraph::V(g)$name <- c(paste0("a", 1:5), paste0("b", 1:4))
  g
}

# A barbell: two cliques joined by a single edge. That edge is the unique bridge
# and the only route between the halves, which makes it the right fixture for
# edge-weight barrier tests and for shortest-path extraction.
barbell_graph <- function(clique_size = 4) {
  g <- igraph::disjoint_union(igraph::make_full_graph(clique_size),
                              igraph::make_full_graph(clique_size))
  igraph::V(g)$name <- c(paste0("L", seq_len(clique_size)),
                         paste0("R", seq_len(clique_size)))
  g <- igraph::add_edges(g, c("L1", "R1"))
  g
}

# robustness_analysis() drives a txtProgressBar when repeating random removals,
# which would otherwise flood the test output.
without_progress <- function(expr) {
  res <- NULL
  utils::capture.output(res <- force(expr))
  res
}

# plot_Net() and highlight_nodes() draw with base graphics. Without a device, R
# silently writes an Rplots.pdf into the test directory; this routes drawing to
# the null device and always closes it.
draw_quietly <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}
