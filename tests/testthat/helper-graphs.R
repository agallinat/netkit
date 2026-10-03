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
