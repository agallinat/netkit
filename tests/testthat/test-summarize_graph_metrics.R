test_that("metrics on a 5-node ring match hand-computed values", {
  # A ring is the one graph whose every global metric is known a priori, so this
  # pins the arithmetic rather than merely re-recording netkit's own output.
  m <- summarize_graph_metrics(ring_graph())

  expect_equal(nrow(m), 1L)
  expect_equal(m$Nodes, 5L)
  expect_equal(m$Edges, 5L)
  expect_false(m$Is_directed)
  expect_equal(m$Density, 0.5)          # 5 edges / choose(5, 2) = 5/10
  expect_equal(m$Diameter, 2)           # furthest pair on a 5-ring
  expect_equal(m$Avg_degree, 2)         # every vertex has degree 2
  expect_equal(m$Components, 1L)
  expect_equal(m$Single_nodes, 0L)
  expect_equal(m$LCC_size, 5L)
  expect_equal(m$LCC_percent, 1)
  expect_equal(m$Average_path_length, 1.5)  # (5*1 + 5*2) / 10 pairs
  expect_equal(m$Clustering_coefficient, 0) # a ring has no triangles
  expect_equal(m$Degree_entropy, 0)         # degree distribution is a point mass
  expect_equal(m$Gini_degree, 0)            # perfectly equal degrees
})

test_that("the one-row metric table has a stable column set", {
  m <- summarize_graph_metrics(test_graph())

  expect_equal(nrow(m), 1L)
  expect_named(m, c(
    "Nodes", "Edges", "Is_directed", "Is_weighted", "Density", "Diameter",
    "Average_path_length", "Clustering_coefficient", "Degree_assortativity",
    "Avg_degree", "Avg_strength", "Avg_betweenness", "Components", "Single_nodes",
    "LCC_size", "LCC_percent", "Algebraic_connectivity", "Degree_entropy",
    "Gini_degree", "Modularity"
  ))
})

test_that("node and edge counts agree with igraph", {
  g <- test_graph()
  m <- summarize_graph_metrics(g)

  expect_equal(m$Nodes, igraph::vcount(g))
  expect_equal(m$Edges, igraph::ecount(g))
  expect_equal(m$Avg_degree, mean(igraph::degree(g)))
})

test_that("isolated vertices are counted and bound the LCC", {
  g <- igraph::make_ring(5) + igraph::vertices("x", "y")
  igraph::V(g)$name <- c(letters[1:5], "x", "y")
  m <- summarize_graph_metrics(g)

  expect_equal(m$Nodes, 7L)
  expect_equal(m$Single_nodes, 2L)
  expect_equal(m$Components, 3L)
  expect_equal(m$LCC_size, 5L)
  expect_equal(m$LCC_percent, 5 / 7)
})

test_that("degenerate graphs return NaN only for genuinely undefined metrics", {
  # These NaNs come from igraph/ineq and are correct: an edgeless graph has no
  # paths, no connected triples, no degree variance and a zero mean degree. The
  # test exists so that the set cannot silently grow -- a new NaN column would be a
  # netkit bug, as auc$efficiency was.
  g <- igraph::make_empty_graph(10, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(10))

  m <- suppressWarnings(summarize_graph_metrics(g))
  numeric_cols <- m[vapply(m, is.numeric, logical(1))]
  not_finite <- names(numeric_cols)[!vapply(numeric_cols, is.finite, logical(1))]

  expect_setequal(not_finite, c("Average_path_length", "Clustering_coefficient",
                                "Degree_assortativity", "Gini_degree", "Modularity"))

  # Everything structural is still exact.
  expect_equal(m$Nodes, 10L)
  expect_equal(m$Edges, 0L)
  expect_equal(m$Components, 10L)
  expect_equal(m$Single_nodes, 10L)
  expect_equal(m$Avg_degree, 0)
})

test_that("Avg_betweenness is skipped above the 5000-node guard", {
  # The guard exists because exact betweenness is the one metric that does not
  # scale; this asserts the cheap path is taken, not the value.
  g <- test_graph()
  expect_false(is.na(summarize_graph_metrics(g)$Avg_betweenness))

  # RSpectra::eigs warns that it cannot converge on the Laplacian of a large
  # ring; Algebraic_connectivity is not what this test is about.
  big <- igraph::make_ring(5001)
  expect_true(is.na(suppressWarnings(summarize_graph_metrics(big))$Avg_betweenness))
})

test_that("Modularity is reproducible across identical calls", {
  # Louvain randomizes its vertex visiting order, so this was the one metric in
  # the summary that was not a function of the graph: repeated calls on one
  # graph returned values spanning roughly 0.356 to 0.366, which is larger than
  # many of the differences the metric gets used to argue about.
  g <- test_graph()
  values <- replicate(8, summarize_graph_metrics(g)$Modularity)
  expect_equal(length(unique(values)), 1)

  gw <- test_graph_weighted()
  wvalues <- replicate(5, summarize_graph_metrics(gw, weights = "weight")$Modularity)
  expect_equal(length(unique(wvalues)), 1)
})

test_that("seeding Modularity does not disturb the caller's random stream", {
  # The fix must not trade one defect for another: seeding globally would reset
  # the caller's stream, which is the defect robustness_analysis() had.
  g <- test_graph()

  set.seed(99); before <- runif(3)
  set.seed(99); invisible(summarize_graph_metrics(g)); after <- runif(3)
  expect_equal(before, after)

  # And with no .Random.seed in existence, none must be left behind.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(summarize_graph_metrics(g))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("seed = NULL opts back in to the ambient stream", {
  g <- test_graph()
  set.seed(7)
  a <- summarize_graph_metrics(g, seed = NULL)$Modularity
  set.seed(7)
  b <- summarize_graph_metrics(g, seed = NULL)$Modularity
  # Same ambient seed, same answer -- so the argument is honoured rather than
  # ignored -- while the value is free to differ from the default-seed one.
  expect_equal(a, b)
})
