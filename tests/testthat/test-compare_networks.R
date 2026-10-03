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

# --- Regression: the three similarity statistics are three statistics -------
#
# `edge_overlap` used to be the identical expression to `jaccard_similarity` --
# the same number reported twice under two names, which looks like corroboration
# and is not. It is now the overlap coefficient, dividing by the smaller edge
# set rather than by the union.

test_that("jaccard_similarity and edge_overlap are no longer the same number", {
  # Graphs of very different size: this is exactly the case where Jaccard and the
  # overlap coefficient diverge, and where reporting one twice is most misleading.
  big <- test_graph_gnp(n = 60, p = 0.12, seed = 11)
  small <- igraph::induced_subgraph(big, igraph::V(big)[1:12])

  sim <- compare_networks(big, small, show_PL = FALSE)$similarity

  expect_false(isTRUE(all.equal(sim$jaccard_similarity, sim$edge_overlap)))
  expect_gt(sim$edge_overlap, sim$jaccard_similarity)
})

test_that("edge_overlap is 1 for a subgraph and Jaccard is not", {
  big <- test_graph_gnp(n = 50, p = 0.1, seed = 12)
  sub <- igraph::induced_subgraph(big, igraph::V(big)[1:20])

  sim <- compare_networks(big, sub, show_PL = FALSE)$similarity

  # Every edge of the subgraph is an edge of the parent, so the smaller set is
  # wholly contained: the overlap coefficient is exactly 1.
  expect_equal(sim$edge_overlap, 1)
  expect_lt(sim$jaccard_similarity, 1)
})

test_that("identical graphs score 1 on all three statistics", {
  g <- test_graph_gnp()
  sim <- compare_networks(g, g, show_PL = FALSE)$similarity
  expect_equal(sim$jaccard_similarity, 1)
  expect_equal(sim$node_overlap, 1)
  expect_equal(sim$edge_overlap, 1)
})

test_that("an empty edge set yields NaN rather than a division error", {
  empty <- igraph::make_empty_graph(n = 5, directed = FALSE)
  igraph::V(empty)$name <- paste0("e", 1:5)
  sim <- suppressWarnings(compare_networks(empty, empty, show_PL = FALSE)$similarity)
  expect_true(is.nan(sim$jaccard_similarity))
  expect_true(is.nan(sim$edge_overlap))
})

test_that("the single-node notice is a suppressible message, not cat() output", {
  g1 <- igraph::add_vertices(test_graph_gnp(n = 30, seed = 13), 3,
                             name = c("iso1", "iso2", "iso3"))
  g2 <- test_graph_gnp(n = 30, seed = 14)

  expect_message(compare_networks(g1, g2, remove_singles = TRUE, show_PL = FALSE),
                 "Single nodes excluded")
  # cat() writes to stdout and cannot be suppressed; message() can.
  expect_silent(suppressMessages(
    compare_networks(g1, g2, remove_singles = TRUE, show_PL = FALSE)
  ))
})
