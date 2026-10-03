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

test_that("non-graph input is rejected with a clear message", {
  for (fn in list(summarize_graph_metrics, plot_CCDF, find_modules,
                  find_hubs, find_bottlenecks, prepare_diffusion)) {
    expect_error(fn("not a graph"), "igraph object or a data.frame")
  }
  expect_error(network_diffusion("not a graph", seed_nodes = "a"),
               "igraph object or a data.frame")
  expect_error(greedy_seed_selection("not a graph", target_nodes = "a"),
               "igraph object or edge list")
})

test_that("calculate_roles() and layout_horizontal_tree() require an igraph object", {
  # Documents a real inconsistency: these two exports do not implement the
  # data.frame branch of the input preamble. If that is ever unified, these
  # expectations should flip to expect_type()/expect_true().
  expect_error(calculate_roles(el), "Must provide a graph object")
  expect_error(layout_horizontal_tree(el), "Must provide a graph object")
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
