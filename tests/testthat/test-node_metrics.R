# node_metrics() is the node-level counterpart to summarize_graph_metrics().
# Where possible these tests assert a property (a bound, a conservation law, an
# identity) or a value known by hand on the 5-node ring, rather than recording
# whatever the code currently emits.

g <- test_graph()

# --- Shape and the shared return contract -----------------------------------

test_that("node_metrics() returns the full plot/result/graph/method shape", {
  res <- node_metrics(g, plot = FALSE)

  expect_type(res, "list")
  expect_true(all(c("plot", "result", "graph", "method") %in% names(res)))
  expect_null(res$plot)
  expect_s3_class(res$result, "tbl_df")
  expect_s3_class(res$graph, "igraph")
  expect_type(res$method, "character")
  expect_length(res$method, 1)
})

test_that("the result covers every vertex exactly once, in graph order", {
  res <- node_metrics(g, plot = FALSE)
  expect_equal(nrow(res$result), igraph::vcount(g))
  expect_false(any(duplicated(res$result$node)))
  # Graph order, not sorted order: the vertex attributes are assigned
  # positionally, so the table must line up with igraph::V(graph).
  expect_identical(res$result$node, igraph::V(g)$name)
})

test_that("the column set follows `metrics`, in the order requested", {
  res <- node_metrics(g, metrics = c("pagerank", "degree", "coreness"),
                      plot = FALSE)
  expect_named(res$result, c("node", "pagerank", "degree", "coreness"))
})

test_that("an unknown metric is rejected eagerly", {
  expect_error(node_metrics(g, metrics = "importance", plot = FALSE))
})

# --- Values known by hand on the ring ---------------------------------------

test_that("metrics match hand-computed values on a 5-node ring", {
  r <- ring_graph()
  res <- node_metrics(r, metrics = c("degree", "strength", "betweenness",
                                     "harmonic", "eigenvector", "pagerank",
                                     "coreness", "clustering", "eccentricity"),
                      normalized = FALSE, plot = FALSE)
  tb <- res$result

  # Every node of a ring is identical, so every metric is constant.
  expect_true(all(tb$degree == 2))
  expect_true(all(tb$strength == 2))
  expect_true(all(tb$coreness == 2))

  # On a 5-ring each node sits on exactly one shortest path between its two
  # 2-hop neighbors, in one direction: unnormalized betweenness is 1.
  expect_true(all(tb$betweenness == 1))

  # Distances from any node are 1, 1, 2, 2, so harmonic centrality is
  # 1 + 1 + 1/2 + 1/2 = 3 and eccentricity is 2.
  expect_true(all(abs(tb$harmonic - 3) < 1e-10))
  expect_true(all(tb$eccentricity == 2))

  # A ring is triangle-free, so local clustering is 0 everywhere.
  expect_true(all(tb$clustering == 0))

  # PageRank is uniform on a vertex-transitive graph and sums to 1.
  expect_true(all(abs(tb$pagerank - 1 / 5) < 1e-8))

  # igraph scales eigenvector centrality so the maximum is 1.
  expect_true(all(abs(tb$eigenvector - 1) < 1e-8))
})

# --- Properties that must hold on any graph ---------------------------------

test_that("degree agrees with the adjacency matrix row sums", {
  res <- node_metrics(g, metrics = "degree", normalized = FALSE, plot = FALSE)
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  expect_equal(res$result$degree, as.numeric(Matrix::rowSums(A)))
})

test_that("coreness never exceeds degree and pagerank sums to one", {
  res <- node_metrics(g, metrics = c("degree", "coreness", "pagerank"),
                      normalized = FALSE, plot = FALSE)
  expect_true(all(res$result$coreness <= res$result$degree))
  expect_equal(sum(res$result$pagerank), 1, tolerance = 1e-8)
})

test_that("normalization keeps the bounded metrics within [0, 1]", {
  res <- node_metrics(g, metrics = c("degree", "betweenness", "harmonic"),
                      normalized = TRUE, plot = FALSE)
  for (m in c("degree", "betweenness", "harmonic")) {
    v <- res$result[[m]]
    v <- v[is.finite(v)]
    expect_true(all(v >= 0 & v <= 1), info = m)
  }
})

test_that("normalization is a positive rescaling, so it preserves the ranking", {
  raw <- node_metrics(g, metrics = "betweenness", normalized = FALSE,
                      plot = FALSE)$result
  nrm <- node_metrics(g, metrics = "betweenness", normalized = TRUE,
                      plot = FALSE)$result
  expect_equal(rank(raw$betweenness), rank(nrm$betweenness))
})

test_that("local clustering is NaN exactly for nodes with fewer than two neighbors", {
  # Not a defect: such a node has no pair of neighbors, so the quantity is
  # undefined rather than zero. Pinned so it is not "fixed" into a 0.
  res <- node_metrics(g, metrics = c("degree", "clustering"),
                      normalized = FALSE, plot = FALSE)
  expect_equal(is.nan(res$result$clustering), res$result$degree < 2)
})

# --- Chaining ---------------------------------------------------------------

test_that("every metric is attached to the returned graph under its own name", {
  metrics <- c("degree", "betweenness", "pagerank", "coreness")
  res <- node_metrics(g, metrics = metrics, plot = FALSE)

  expect_true(all(metrics %in% igraph::vertex_attr_names(res$graph)))
  for (m in metrics) {
    expect_equal(igraph::vertex_attr(res$graph, m), res$result[[m]], info = m)
  }
  # Pre-existing attributes survive, so the graph keeps chaining.
  expect_true(all(c("name", "score") %in% igraph::vertex_attr_names(res$graph)))
})

test_that("the annotated graph feeds robustness_analysis() directly", {
  # robustness_analysis() already documents "any numeric vertex attribute" as a
  # removal strategy; node_metrics() is what makes that convenient.
  res <- node_metrics(g, metrics = c("pagerank", "coreness"), plot = FALSE)
  rob <- robustness_analysis(res$graph, removal_strategy = "pagerank",
                             steps = 8, metrics = "lcc_size", plot = FALSE)
  expect_true(is.numeric(rob$auc$lcc_size))
})

# --- Weights ----------------------------------------------------------------

test_that("unit weights reproduce the unweighted metrics", {
  gw <- test_graph_unit_weights()
  gu <- igraph::delete_edge_attr(gw, "weight")
  metrics <- c("degree", "strength", "betweenness", "harmonic", "eigenvector",
               "pagerank", "coreness")

  a <- node_metrics(gw, metrics = metrics, weights = "weight", plot = FALSE)$result
  b <- node_metrics(gu, metrics = metrics, plot = FALSE)$result
  expect_equal(a, b, tolerance = 1e-8)
})

test_that("strength equals degree when unweighted and differs when weighted", {
  res_u <- node_metrics(g, metrics = c("degree", "strength"),
                        normalized = FALSE, plot = FALSE)$result
  expect_equal(res_u$strength, res_u$degree)

  gw <- test_graph_weighted()
  res_w <- node_metrics(gw, metrics = c("degree", "strength"),
                        normalized = FALSE, weights = "weight",
                        plot = FALSE)$result
  expect_false(isTRUE(all.equal(res_w$strength, res_w$degree)))
  # Strengths are sums of weights below 1, so they fall under the counts here.
  expect_true(all(res_w$strength <= res_w$degree))
})

# --- The costly-metric guard ------------------------------------------------

test_that("metrics above max_nodes are skipped loudly and returned as NA", {
  # The column set must not depend on the graph's size, or a script written
  # against a small graph breaks silently on a large one.
  expect_warning(
    res <- node_metrics(g, metrics = c("degree", "betweenness"),
                        max_nodes = 10, plot = FALSE),
    "were skipped|was skipped"
  )
  expect_named(res$result, c("node", "degree", "betweenness"))
  expect_true(all(is.na(res$result$betweenness)))
  expect_false(any(is.na(res$result$degree)))
  expect_match(res$method, "skipped above max_nodes")
})

test_that("cheap metrics are unaffected by max_nodes", {
  expect_no_warning(
    node_metrics(g, metrics = c("degree", "coreness", "pagerank"),
                 max_nodes = 10, plot = FALSE)
  )
})

# --- Disconnected graphs ----------------------------------------------------

test_that("the defaults are silent on a disconnected graph", {
  expect_no_warning(node_metrics(test_graph_disconnected(), plot = FALSE))
})

test_that("closeness ranks the smaller component higher, which is why it is not a default", {
  # igraph averages distance over *reachable* vertices only, so on a fragmented
  # graph closeness rewards being in a small component: everything nearby.
  # test_graph_disconnected() is a 5-clique plus a 4-clique, so the 4-clique
  # nodes have 3 neighbors at distance 1 (closeness 1/3) and the 5-clique nodes
  # have 4 (closeness 1/4). Ranking by closeness puts the periphery on top.
  g <- test_graph_disconnected()
  res <- node_metrics(g, metrics = c("closeness", "harmonic"), plot = FALSE,
                      normalized = FALSE)
  big <- res$result$node %in% paste0("a", 1:5)

  # Comparing the extremes, not the vectors: the two components have different
  # sizes, so an elementwise comparison would recycle.
  expect_gt(min(res$result$closeness[!big]), max(res$result$closeness[big]))

  # harmonic centrality does not have that failure mode: unreachable pairs
  # contribute zero rather than being dropped, so the larger component wins.
  expect_gt(min(res$result$harmonic[big]), max(res$result$harmonic[!big]))
  expect_true(all(is.finite(res$result$harmonic)))
})

# --- Plots ------------------------------------------------------------------

test_that("both plot types build without drawing", {
  corr <- node_metrics(g, plot = TRUE, plot_type = "correlation")$plot
  rank <- node_metrics(g, plot = TRUE, plot_type = "ranking", top_n = 5)$plot

  expect_s3_class(corr, "ggplot")
  expect_s3_class(rank, "ggplot")
  # ggplot is lazy, so building is what actually exercises the layers.
  expect_no_error(ggplot2::ggplot_build(corr))
  expect_no_error(ggplot2::ggplot_build(rank))
})

test_that("the correlation heatmap is square, symmetric and diagonal-1", {
  metrics <- c("degree", "betweenness", "pagerank")
  p <- node_metrics(g, metrics = metrics, plot = TRUE,
                    plot_type = "correlation")$plot
  df <- p$data

  expect_equal(nrow(df), length(metrics)^2)
  diag_vals <- df$correlation[as.character(df$metric_x) == as.character(df$metric_y)]
  expect_equal(diag_vals, rep(1, length(metrics)), tolerance = 1e-8)
  expect_true(all(df$correlation >= -1 & df$correlation <= 1))
})

test_that("the correlation plot is NULL when fewer than two metrics vary", {
  # A ring has constant everything, so there is nothing to correlate. Returning
  # NULL beats emitting a cor() warning the caller cannot act on.
  r <- ring_graph()
  expect_no_warning(
    p <- node_metrics(r, metrics = c("degree", "coreness"), plot = TRUE)$plot
  )
  expect_null(p)
})

test_that("the ranking plot shows at most top_n bars per metric", {
  metrics <- c("degree", "pagerank")
  p <- node_metrics(g, metrics = metrics, plot = TRUE,
                    plot_type = "ranking", top_n = 4)$plot
  counts <- table(p$data$metric)
  expect_true(all(counts <= 4))
  expect_setequal(names(counts), metrics)
})

# --- Input contract ---------------------------------------------------------

test_that("node_metrics() honours the shared graph-input contract", {
  expect_error(node_metrics("not a graph"),
               "Input 'graph' must be either an igraph object or a data.frame")
  expect_type(node_metrics(as_edge_df(g), plot = FALSE), "list")
})
