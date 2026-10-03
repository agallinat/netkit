# netkit's analysis functions share a return vocabulary -- `plot`, `result`,
# `graph`, `method` -- and that shared shape is what lets them chain. These
# tests pin the contract so a refactor cannot quietly change it.

g <- test_graph()

test_that("node-classification functions return the full plot/result/graph/method shape", {
  for (res in list(find_hubs(g, plot = FALSE), find_bottlenecks(g, plot = FALSE))) {
    expect_type(res, "list")
    expect_true(all(c("plot", "result", "graph", "method") %in% names(res)))
    expect_s3_class(res$result, "data.frame")
    expect_s3_class(res$graph, "igraph")
    expect_type(res$method, "character")
    expect_length(res$method, 1)
  }
})

test_that("`plot = FALSE` suppresses the plot object", {
  expect_null(find_hubs(g, plot = FALSE)$plot)
  expect_null(find_bottlenecks(g, plot = FALSE)$plot)
})

test_that("`plot = TRUE` returns a plot object rather than drawing one", {
  # find_hubs()/find_bottlenecks() pass their scatterplot through
  # ggExtra::ggMarginal() to add marginal distributions, which returns an
  # already-assembled gtable rather than a ggplot. Assembling that gtable needs a
  # graphics device, so these go through draw_quietly() even though nothing is
  # being plotted on purpose.
  expect_s3_class(draw_quietly(find_hubs(g, plot = TRUE))$plot, "ggExtraPlot")
  expect_s3_class(draw_quietly(find_bottlenecks(g, plot = TRUE))$plot, "ggExtraPlot")

  expect_s3_class(calculate_roles(g, cluster.method = "louvain", plot = TRUE)$plot, "ggplot")
})

test_that("find_modules() has no plot slot because it draws via plot_Net()", {
  # find_modules() delegates to plot_Net(), which uses base plot.igraph, so there
  # is no ggplot object to hand back under either setting.
  expect_false("plot" %in% names(find_modules(g, plot = FALSE)))
  expect_false("plot" %in% names(draw_quietly(find_modules(g, plot = TRUE))))
})

test_that("calculate_roles() and robustness_analysis() omit `plot` entirely when plot = FALSE", {
  # Note the inconsistency with find_hubs()/find_bottlenecks(), which keep the
  # name and set it to NULL. Recorded rather than endorsed: if these are
  # harmonised to always return a `plot` slot, update these expectations.
  roles <- calculate_roles(g, cluster.method = "louvain", plot = FALSE)
  expect_false("plot" %in% names(roles))

  rob <- robustness_analysis(g, removal_strategy = "degree", steps = 5,
                             n_reps = 1, plot = FALSE, seed = 1)
  expect_false("plot" %in% names(rob))
})

test_that("the annotated graph carries the new vertex attribute and chains onward", {
  hubs <- find_hubs(g, plot = FALSE)
  expect_true("is_hub" %in% igraph::vertex_attr_names(hubs$graph))
  expect_type(igraph::V(hubs$graph)$is_hub, "logical")

  bottle <- find_bottlenecks(hubs$graph, plot = FALSE)
  expect_true(all(c("is_hub", "is_bottleneck") %in%
                    igraph::vertex_attr_names(bottle$graph)))

  mods <- find_modules(bottle$graph, plot = FALSE)
  expect_true(all(c("is_hub", "is_bottleneck", "module") %in%
                    igraph::vertex_attr_names(mods$graph)))

  # Pre-existing attributes survive the round trip.
  expect_true("score" %in% igraph::vertex_attr_names(mods$graph))
})

test_that("node-level result tables cover the graph's vertices exactly once", {
  for (res in list(find_hubs(g, plot = FALSE), find_bottlenecks(g, plot = FALSE))) {
    expect_setequal(res$result$node, igraph::V(g)$name)
    expect_false(any(duplicated(res$result$node)))
  }
})
