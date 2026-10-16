g1 <- test_graph()
g2 <- test_graph_gnp()

test_that("compare_networks() returns the documented structure", {
  res <- compare_networks(g1, g2)

  expect_named(res, c("plot", "global_topology", "similarity", "edge_enrichment",
                      "ks_test"))
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

# --- Shared edges against chance --------------------------------------------
#
# A shared-edge count alone says nothing: two dense networks over one vertex set
# overlap by arithmetic. These assert the arithmetic of the expectation from
# first principles rather than re-recording what the function emits.

test_that("edge_enrichment reports the documented columns", {
  ee <- compare_networks(g1, g2, show_PL = FALSE)$edge_enrichment

  expect_s3_class(ee, "data.frame")
  expect_equal(nrow(ee), 1L)
  expect_named(ee, c("n_shared_nodes", "n_possible_edges", "edges1", "edges2",
                     "shared_edges", "expected_shared_edges",
                     "fold_enrichment", "p_value"))
})

test_that("the expectation is m1 * m2 / choose(n, 2) over the shared vertices", {
  # g1 and g2 share no vertex names, so the pool would be empty. Two graphs over
  # one vertex set are what the statistic is for.
  a <- test_graph_gnp(n = 50, p = 0.08, seed = 21)
  b <- test_graph_gnp(n = 50, p = 0.08, seed = 22)
  ee <- compare_networks(a, b, show_PL = FALSE)$edge_enrichment

  # Recomputed here from the reported draw sizes and pool, which is the identity
  # the column claims to be. Deriving it independently of the function is the
  # point: the defect this guards against is a count and a denominator coming
  # from two different universes.
  expect_equal(ee$expected_shared_edges,
               ee$edges1 * ee$edges2 / ee$n_possible_edges)
  expect_equal(ee$n_possible_edges, choose(ee$n_shared_nodes, 2))
  expect_equal(ee$fold_enrichment, ee$shared_edges / ee$expected_shared_edges)
  expect_true(ee$p_value >= 0 && ee$p_value <= 1)
})

test_that("the draw sizes count only edges inside the shared vertex set", {
  # The halves share vertices s1..s6 and nothing else. Every edge touching a
  # private vertex is outside the pool and must not enter m1 or m2, even though
  # ecount() counts it.
  gA <- igraph::graph_from_edgelist(
    cbind(c("s1", "s2", "s3", "a1", "a2"), c("s2", "s3", "s4", "s1", "a1")),
    directed = FALSE)
  gB <- igraph::graph_from_edgelist(
    cbind(c("s1", "s2", "s5", "b1", "s4"), c("s2", "s3", "s6", "s1", "b1")),
    directed = FALSE)

  ee <- compare_networks(gA, gB, show_PL = FALSE)$edge_enrichment

  # A holds s1..s4, B holds s1..s6 minus nothing relevant: the shared set is
  # s1..s4. a1, a2, s5, s6 and b1 are each in one network only.
  expect_equal(ee$n_shared_nodes, 4L)
  expect_equal(ee$n_possible_edges, choose(4, 2))
  # A's qualifying edges: s1-s2, s2-s3, s3-s4. The a1/a2 edges are excluded.
  expect_equal(ee$edges1, 3L)
  expect_lt(ee$edges1, igraph::ecount(gA))
  # B's: s1-s2, s2-s3 only; s5-s6, b1-s1 and s4-b1 leave the shared set.
  expect_equal(ee$edges2, 2L)
  expect_lt(ee$edges2, igraph::ecount(gB))
  expect_equal(ee$shared_edges, 2L)
})

test_that("a graph against itself is strongly enriched; unrelated graphs are not", {
  same <- compare_networks(g1, g1, show_PL = FALSE)$edge_enrichment
  expect_equal(same$shared_edges, same$edges1)
  expect_gt(same$fold_enrichment, 1)
  expect_lt(same$p_value, 0.01)

  # Two independent draws over the same vertices: the observed overlap should sit
  # near its expectation, so the test must not call that enrichment. Asserted as
  # a band rather than a value, since the graphs are random.
  n <- 120
  ind1 <- igraph::sample_gnp(n, 0.06); igraph::V(ind1)$name <- paste0("v", 1:n)
  ind2 <- igraph::sample_gnp(n, 0.06); igraph::V(ind2)$name <- paste0("v", 1:n)
  ee <- compare_networks(ind1, ind2, show_PL = FALSE)$edge_enrichment

  expect_gt(ee$fold_enrichment, 0.4)
  expect_lt(ee$fold_enrichment, 2.5)
  expect_gt(ee$p_value, 0.01)
})

test_that("a nested subgraph is enriched well above chance", {
  big <- test_graph_gnp(n = 60, p = 0.12, seed = 11)
  sub <- igraph::induced_subgraph(big, igraph::V(big)[1:20])

  ee <- compare_networks(big, sub, show_PL = FALSE)$edge_enrichment

  # Every edge of the subgraph is an edge of the parent, restricted to the same
  # 20 vertices, so the overlap is total and the pool is choose(20, 2).
  expect_equal(ee$n_shared_nodes, 20L)
  expect_equal(ee$shared_edges, ee$edges2)
  expect_gt(ee$fold_enrichment, 1)
  expect_lt(ee$p_value, 0.05)
})

test_that("direction, multiplicity and self-loops are collapsed before counting", {
  # a-a is a loop, b-a repeats a-b: neither is a pair in choose(n, 2), so the
  # draw size is 2 (a|b and c|d) rather than ecount()'s 4.
  m1 <- igraph::graph_from_edgelist(
    cbind(c("a", "a", "b", "c"), c("a", "b", "a", "d")), directed = FALSE)
  m2 <- igraph::graph_from_edgelist(cbind(c("a", "c"), c("b", "d")),
                                    directed = FALSE)

  ee <- compare_networks(m1, m2, show_PL = FALSE)$edge_enrichment

  expect_equal(ee$edges1, 2L)
  expect_equal(ee$edges2, 2L)
  expect_equal(ee$shared_edges, 2L)
  expect_equal(ee$n_possible_edges, choose(4, 2))
  expect_equal(ee$expected_shared_edges, 2 * 2 / choose(4, 2))
})

test_that("an empty pool gives NaN and NA rather than a zero-denominator number", {
  g <- test_graph_gnp(n = 20, seed = 3)
  disjoint <- test_graph_gnp(n = 20, seed = 4)
  igraph::V(disjoint)$name <- paste0("z", seq_len(igraph::vcount(disjoint)))

  ee <- compare_networks(g, disjoint, show_PL = FALSE)$edge_enrichment

  expect_equal(ee$n_shared_nodes, 0L)
  expect_equal(ee$n_possible_edges, 0)
  expect_true(is.nan(ee$expected_shared_edges))
  expect_true(is.nan(ee$fold_enrichment))
  expect_true(is.na(ee$p_value))
})

test_that("an edgeless pair gives a zero expectation and p = 1", {
  empty <- igraph::make_empty_graph(n = 5, directed = FALSE)
  igraph::V(empty)$name <- paste0("e", 1:5)
  ee <- suppressWarnings(
    compare_networks(empty, empty, show_PL = FALSE))$edge_enrichment

  # The pool exists but neither network drew from it. Nothing was observed and
  # nothing was expected, so there is no enrichment to report -- but the p-value
  # is well defined, and is 1.
  expect_equal(ee$n_possible_edges, choose(5, 2))
  expect_equal(ee$expected_shared_edges, 0)
  expect_true(is.nan(ee$fold_enrichment))
  expect_equal(ee$p_value, 1)
})

test_that("unnamed graphs key edges by index in the same way as before", {
  # as_edgelist() falls back to indices for an unnamed graph and the function now
  # resolves them itself; the set statistics must be unchanged by that.
  a <- igraph::sample_gnp(40, 0.1); b <- igraph::sample_gnp(40, 0.1)
  sim_unnamed <- compare_networks(a, b, show_PL = FALSE)$similarity

  igraph::V(a)$name <- as.character(seq_len(40))
  igraph::V(b)$name <- as.character(seq_len(40))
  sim_named <- compare_networks(a, b, show_PL = FALSE)$similarity

  expect_equal(sim_unnamed$jaccard_similarity, sim_named$jaccard_similarity)
  expect_equal(sim_unnamed$edge_overlap, sim_named$edge_overlap)
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
