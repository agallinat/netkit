# netkit's analysis functions share a return vocabulary -- `plot`, `result`,
# `graph`, `method` -- and that shared shape is what lets them chain. These
# tests pin the contract so a refactor cannot quietly change it.

g <- test_graph()

test_that("node-classification functions return the full plot/result/graph/method shape", {
  # calculate_roles() is included here deliberately: it used to return
  # plot/roles_definitions/result only, with no `graph` and no `method`, and
  # nothing in this file noticed. Listing it by name is what stops that
  # recurring.
  classifiers <- list(
    find_hubs        = find_hubs(g, plot = FALSE),
    find_bottlenecks = find_bottlenecks(g, plot = FALSE),
    calculate_roles  = calculate_roles(g, cluster.method = "louvain", plot = FALSE)
  )
  for (nm in names(classifiers)) {
    res <- classifiers[[nm]]
    expect_type(res, "list")
    # Classed for printing, but still a list in every other respect. See
    # test-print.R.
    expect_s3_class(res, "netkit_result")
    expect_true(all(c("plot", "result", "graph", "method") %in% names(res)), info = nm)
    expect_s3_class(res$result, "data.frame")
    expect_s3_class(res$graph, "igraph")
    expect_type(res$method, "character")
    expect_length(res$method, 1)
  }
})

test_that("every analysis function exposes a `result` table", {
  # `result` is the package-wide name for the node- or step-level table.
  # find_modules() and robustness_analysis() historically used `module_table`
  # and `summary` instead; both now also carry `result`, so the vocabulary has
  # no exceptions left to remember.
  results <- list(
    find_hubs           = find_hubs(g, plot = FALSE),
    find_bottlenecks    = find_bottlenecks(g, plot = FALSE),
    calculate_roles     = calculate_roles(g, cluster.method = "louvain", plot = FALSE),
    find_modules        = find_modules(g, plot = FALSE),
    robustness_analysis = robustness_analysis(g, removal_strategy = "degree",
                                              steps = 5, plot = FALSE)
  )
  for (nm in names(results)) {
    expect_true("result" %in% names(results[[nm]]), info = nm)
    expect_s3_class(results[[nm]]$result, "data.frame")
  }

  # The deprecated aliases still point at the same object.
  expect_identical(results$find_modules$result, results$find_modules$module_table)
  expect_identical(results$robustness_analysis$result,
                   results$robustness_analysis$summary)
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

test_that("every plot-returning function keeps a `plot` element in both modes", {
  # The return shape must not depend on the arguments: `plot` is always present,
  # and NULL when plot = FALSE. calculate_roles() and robustness_analysis() used
  # to drop the name entirely, which made names() and str() vary by call.
  off <- list(
    find_hubs             = find_hubs(g, plot = FALSE),
    find_bottlenecks      = find_bottlenecks(g, plot = FALSE),
    calculate_roles       = calculate_roles(g, cluster.method = "louvain", plot = FALSE),
    robustness_analysis   = robustness_analysis(g, removal_strategy = "degree",
                                                steps = 5, n_reps = 1,
                                                plot = FALSE, seed = 1),
    greedy_seed_selection = greedy_seed_selection(g, target_nodes = c("n10", "n11"),
                                                  k = 2, plot = FALSE)
  )
  for (nm in names(off)) {
    expect_true("plot" %in% names(off[[nm]]), info = nm)
    expect_null(off[[nm]]$plot, info = nm)
  }

  on <- list(
    find_hubs           = draw_quietly(find_hubs(g, plot = TRUE)),
    find_bottlenecks    = draw_quietly(find_bottlenecks(g, plot = TRUE)),
    calculate_roles     = calculate_roles(g, cluster.method = "louvain", plot = TRUE),
    robustness_analysis = robustness_analysis(g, removal_strategy = "degree",
                                              steps = 5, n_reps = 1,
                                              plot = TRUE, seed = 1)
  )
  for (nm in names(on)) {
    expect_true("plot" %in% names(on[[nm]]), info = nm)
    expect_false(is.null(on[[nm]]$plot), info = nm)
  }
})

test_that("the plot element is named `plot` everywhere, including compare_networks()", {
  # compare_networks() used to call it CCDF_plot, which broke the shared
  # vocabulary. Renamed while netkit is still unreleased and has no users.
  res <- compare_networks(g, test_graph_gnp())
  expect_true("plot" %in% names(res))
  expect_s3_class(res$plot, "ggplot")
  expect_false("CCDF_plot" %in% names(res))
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
