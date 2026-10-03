g1 <- test_graph()
g2 <- test_graph_gnp()

test_that("compare_networks() returns the documented structure", {
  res <- compare_networks(g1, g2)

  expect_named(res, c("plot", "global_topology", "similarity", "ks_test"))
  expect_s3_class(res$plot, "ggplot")
  expect_s3_class(res$global_topology, "data.frame")
})

test_that("global_topology stacks one row per input network", {
  res <- compare_networks(g1, g2)

  expect_equal(nrow(res$global_topology), 2L)
  # It is built by rbind()-ing summarize_graph_metrics() output, so the columns
  # must stay in step with that function.
  expect_equal(names(res$global_topology),
               names(summarize_graph_metrics(g1)))
  expect_equal(res$global_topology$Nodes,
               c(igraph::vcount(g1), igraph::vcount(g2)))
  expect_equal(res$global_topology$Edges,
               c(igraph::ecount(g1), igraph::ecount(g2)))
})

test_that("similarity reports the three documented overlap measures", {
  # The @return block used to name elements that did not exist (metrics1,
  # metrics2, jaccard_similarity at top level); these are the real ones.
  res <- compare_networks(g1, g2)

  expect_s3_class(res$similarity, "data.frame")
  expect_equal(nrow(res$similarity), 1L)
  expect_named(res$similarity, c("jaccard_similarity", "node_overlap", "edge_overlap"))
  for (col in names(res$similarity)) {
    expect_true(res$similarity[[col]] >= 0 && res$similarity[[col]] <= 1, info = col)
  }
})

test_that("a graph compared with itself overlaps completely", {
  res <- compare_networks(g1, g1)
  expect_equal(res$similarity$jaccard_similarity, 1)
  expect_equal(res$similarity$node_overlap, 1)
  expect_equal(res$similarity$edge_overlap, 1)
})

test_that("the KS test compares the two degree distributions", {
  res <- compare_networks(g1, g2)

  expect_s3_class(res$ks_test, "htest")
  expect_true(res$ks_test$p.value >= 0 && res$ks_test$p.value <= 1)
})

test_that("comparing a graph with itself reports maximal similarity", {
  res <- compare_networks(g1, g1)

  expect_equal(res$global_topology$Nodes, c(igraph::vcount(g1), igraph::vcount(g1)))
  # Identical degree distributions: KS cannot reject the null.
  expect_equal(res$ks_test$statistic, c(D = 0), ignore_attr = TRUE)
  expect_equal(res$ks_test$p.value, 1)
})

test_that("compare_networks() options are accepted", {
  expect_type(compare_networks(g1, g2, show_PL = FALSE), "list")
  expect_type(compare_networks(g1, g2, remove_singles = TRUE), "list")
})
