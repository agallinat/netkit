# The "flexible graph input" preamble is netkit's most-repeated convention:
# accept an igraph object OR a data.frame edge list, error otherwise, then
# backfill V(graph)$name when absent. It is copy-pasted across ~12 functions, so
# these tests exist to keep every copy in agreement -- and to record the two
# exports that do NOT follow it.

g  <- test_graph()
el <- as_edge_df(g)

test_that("edge-list data.frames are accepted wherever an igraph is", {
  expect_s3_class(summarize_graph_metrics(el), "data.frame")
  expect_s3_class(plot_CCDF(el), "ggplot")
  expect_type(find_modules(el, plot = FALSE), "list")
  expect_type(find_hubs(el, plot = FALSE), "list")
  expect_type(find_bottlenecks(el, plot = FALSE), "list")
  expect_type(prepare_diffusion(el, method = "laplacian"), "list")
  expect_s3_class(network_diffusion(el, seed_nodes = "n1"), "data.frame")
  expect_type(
    robustness_analysis(el, removal_strategy = "degree", steps = 5,
                        n_reps = 1, plot = FALSE, seed = 1),
    "list"
  )
})

test_that("non-graph input is rejected with one consistent message", {
  # All of these route through the shared as_netkit_graph() validator, so the
  # wording is identical rather than varying per function as it used to.
  for (fn in list(summarize_graph_metrics, plot_CCDF, find_modules, find_hubs,
                  find_bottlenecks, prepare_diffusion, plot_Net, assign_attributes,
                  robustness_analysis, highlight_nodes)) {
    expect_error(fn("not a graph"),
                 "Input 'graph' must be either an igraph object or a data.frame")
  }
  expect_error(network_diffusion("not a graph", seed_nodes = "a"),
               "Input 'graph' must be either an igraph object or a data.frame")
  expect_error(greedy_seed_selection("not a graph", target_nodes = "a"),
               "Input 'graph' must be either an igraph object or a data.frame")
})

test_that("the offending argument is named for multi-graph functions", {
  g <- test_graph()
  expect_error(compare_networks("not a graph", g), "Input 'graph1' must be")
  expect_error(compare_networks(g, "not a graph"), "Input 'graph2' must be")
})

test_that("calculate_roles() and layout_horizontal_tree() accept edge lists too", {
  # These two were the last exports not implementing the data.frame branch of the
  # input preamble; both now route through as_netkit_graph() like the rest.
  roles <- calculate_roles(el, cluster.method = "louvain", plot = FALSE)
  expect_type(roles, "list")
  expect_setequal(roles$result$node, igraph::V(g)$name)

  lay <- layout_horizontal_tree(el)
  expect_true(is.matrix(lay))
  expect_equal(nrow(lay), igraph::vcount(g))

  # And they reject non-graph input with the same shared message as everything else.
  expect_error(calculate_roles("not a graph"),
               "Input 'graph' must be either an igraph object or a data.frame")
  expect_error(layout_horizontal_tree("not a graph"),
               "Input 'graph' must be either an igraph object or a data.frame")
})

test_that("missing vertex names are backfilled with indices", {
  unnamed <- igraph::sample_gnp(20, 0.15, directed = FALSE)
  expect_null(igraph::V(unnamed)$name)

  res <- network_diffusion(unnamed, seed_nodes = "1", method = "rwr")
  expect_setequal(res$node, as.character(seq_len(20)))

  # The input graph is not mutated in the caller's frame.
  expect_null(igraph::V(unnamed)$name)
})

test_that("directed graphs are accepted and coerced where documented", {
  dg <- test_graph_directed()
  expect_true(igraph::is_directed(dg))

  metrics <- summarize_graph_metrics(dg)
  expect_true(metrics$Is_directed)
  # Is_directed records the *input*, while the metrics themselves are computed
  # on the collapsed undirected graph.
  expect_false(metrics$Edges > igraph::ecount(dg))

  expect_type(robustness_analysis(dg, removal_strategy = "degree", steps = 5,
                                  n_reps = 1, plot = FALSE, seed = 1), "list")
})
