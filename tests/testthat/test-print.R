# netkit's analysis functions return lists holding an igraph and a full result
# table. Printing one with R's default method dumps all of it: a 100-graph
# null_model() ensemble ran to 1404 lines at the console. These tests pin the
# summarized printing, and -- more importantly -- pin that classing the returns
# did not change their list semantics.

g <- test_graph(n = 40)

test_that("printing a result is a short summary, not a dump of the graph", {
  res <- without_empty_selection(find_hubs(g, plot = FALSE))
  out <- capture.output(print(res))

  # The default method prints the whole 40-row tibble and the igraph edge list.
  # Anything near that length means the method is not being dispatched.
  expect_lt(length(out), 20)
  expect_match(out[1], "netkit result: find_hubs()", fixed = TRUE)

  # The three things worth seeing without asking: what was computed, the head
  # of the table, and the shape of the annotated graph.
  expect_true(any(grepl("^Method: Hub nodes identified", out)))
  expect_true(any(grepl("$result", out, fixed = TRUE)))
  expect_true(any(grepl("40 nodes", out, fixed = TRUE)))
})

test_that("a null ensemble prints its size rather than its graphs", {
  # The regression this guards: print(null_model(g)) used to emit 14 lines per
  # graph, so the documented default of n = 100 produced over 1400.
  nulls <- null_model(g, n = 20, seed = 1)
  out <- capture.output(print(nulls))

  expect_lt(length(out), 6)
  expect_match(out[1], "20 graphs, model 'rewire'", fixed = TRUE)
  expect_true(any(grepl("40 nodes", out, fixed = TRUE)))
})

test_that("subsetting a null ensemble keeps it an ensemble", {
  nulls <- null_model(g, n = 6, seed = 1)
  sub <- nulls[1:3]

  expect_s3_class(sub, "netkit_null")
  expect_length(sub, 3)
  expect_identical(attr(sub, "model"), attr(nulls, "model"))
  expect_lt(length(capture.output(print(sub))), 6)

  # And it is still usable as a null ensemble.
  expect_identical(sub[[1]], nulls[[1]])
})

test_that("a diffusion kernel prints its shape rather than its matrices", {
  # The kernel holds a Laplacian, a Cholesky factor and a transition matrix.
  kernel <- prepare_diffusion(g, method = "rwr")
  out <- capture.output(print(kernel))

  expect_lt(length(out), 6)
  expect_match(out[1], "method 'rwr'", fixed = TRUE)
  expect_true(any(grepl("40 x 40", out, fixed = TRUE)))

  # Classing it must not stop network_diffusion() from reusing it.
  scores <- network_diffusion(g, seed_nodes = c("n1", "n2"), method = "rwr",
                              precompute = kernel)
  expect_s3_class(scores, "data.frame")
})

test_that("classing the returns leaves their list semantics untouched", {
  res <- without_empty_selection(find_hubs(g, plot = FALSE))

  expect_type(res, "list")
  expect_true(is.list(res))
  expect_named(res, c("plot", "method", "result", "graph"))
  expect_s3_class(res$result, "data.frame")
  expect_s3_class(res$graph, "igraph")
  expect_null(res$plot)

  # unclass() is the documented escape hatch back to default printing.
  expect_identical(class(unclass(res)), "list")
  expect_identical(unclass(res)$result, res$result)
})

test_that("print returns its argument invisibly", {
  res <- without_empty_selection(find_hubs(g, plot = FALSE))
  capture.output(expect_invisible(print(res)))

  capture.output(out <- print(res))
  expect_identical(out, res)

  nulls <- null_model(g, n = 2, seed = 1)
  capture.output(out_null <- print(nulls))
  expect_identical(out_null, nulls)
})

test_that("results with no table still print their elements", {
  # compare_networks() carries neither `method` nor `result`, so the printer has
  # nothing to lead with. It must still list what is there rather than printing
  # an empty header.
  cmp <- capture.output(print(compare_networks(g, test_graph_gnp())))
  expect_lt(length(cmp), 12)
  expect_true(any(grepl("$similarity", cmp, fixed = TRUE)))
  expect_true(any(grepl("$edge_enrichment", cmp, fixed = TRUE)))

  # greedy_seed_selection() used to be the other example here. It now carries
  # both a `method` and a `result` table, so it belongs with the normal case
  # below; compare_networks() is the only remaining return with neither.
})

test_that("greedy_seed_selection() prints its method, its table and its vectors", {
  res <- without_leaf_warning(
    greedy_seed_selection(g, target_nodes = c("n10", "n11"), k = 2,
                          plot = FALSE))
  out <- capture.output(print(res))

  # Three sections, and still a summary rather than a dump: the default method
  # would print the whole tibble plus four vectors with no structure.
  expect_lt(length(out), 25)
  expect_true(any(grepl("^Method: Greedy seed selection", out)))
  expect_true(any(grepl("$result", out, fixed = TRUE)))
  expect_true(any(grepl("$selected_seeds", out, fixed = TRUE)))
  expect_true(any(grepl("$final_target_score", out, fixed = TRUE)))
  # The per-seed diagnostic columns are what the table is for, so they have to
  # survive the summarized printing.
  expect_true(any(grepl("gain_z", out, fixed = TRUE)))
})

test_that("a present plot is described, not drawn", {
  # The printer must never trigger a device: describing the slot is the whole
  # point. find_hubs() carries a ggExtra gtable rather than a ggplot, so both
  # spellings have to be recognized.
  res <- without_empty_selection(draw_quietly(find_hubs(g, plot = TRUE)))
  out <- capture.output(print(res))

  expect_true(any(grepl("print(x$plot) to draw", out, fixed = TRUE)))
  expect_lt(length(out), 20)
})

test_that("every list-returning function carries the netkit_result class", {
  # Listing them by name is what stops a new function from quietly skipping
  # as_netkit_result() and printing as a raw list.
  results <- without_empty_selection(list(
    find_hubs             = find_hubs(g, plot = FALSE),
    find_bottlenecks      = find_bottlenecks(g, plot = FALSE),
    calculate_roles       = calculate_roles(g, cluster.method = "louvain",
                                            plot = FALSE),
    find_modules          = find_modules(g, plot = FALSE),
    node_metrics          = node_metrics(g, plot = FALSE),
    compare_networks      = compare_networks(g, test_graph_gnp()),
    extract_subnetwork    = extract_subnetwork(g, c("n1", "n20"),
                                               method = "shortest_paths",
                                               plot = FALSE),
    metric_significance   = metric_significance(g, metrics = "Modularity",
                                                n = 5, seed = 1, plot = FALSE),
    small_worldness       = small_worldness(g, n = 5, seed = 1),
    robustness_analysis   = robustness_analysis(g, removal_strategy = "degree",
                                                steps = 5, plot = FALSE),
    greedy_seed_selection = greedy_seed_selection(g, target_nodes = "n10",
                                                  k = 1, plot = FALSE)
  ))

  for (nm in names(results)) {
    expect_s3_class(results[[nm]], "netkit_result")
    # Each result also carries a per-function subclass, so a caller can
    # dispatch on one kind of result without dispatching on all of them.
    expect_s3_class(results[[nm]], paste0("netkit_", nm))
    expect_identical(attr(results[[nm]], "netkit_fn"), nm)
    expect_lt(length(capture.output(print(results[[nm]]))), 25)
  }
})
