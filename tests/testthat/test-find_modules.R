g <- test_graph()

test_that("find_modules() returns the documented structure", {
  set.seed(1)
  res <- find_modules(g, plot = FALSE)

  expect_named(res, c("result", "module_table", "n_modules", "modularity",
                      "subgraphs", "method", "graph"))
  expect_s3_class(res$module_table, "tbl_df")
  expect_named(res$module_table, c("node", "module"))
  expect_type(res$n_modules, "integer")
  expect_gt(res$n_modules, 0)
  expect_equal(res$method, "louvain")
  expect_s3_class(res$graph, "igraph")
  expect_true("module" %in% igraph::vertex_attr_names(res$graph))
})

test_that("module assignments are consistent with n_modules", {
  set.seed(1)
  res <- find_modules(g, plot = FALSE)

  expect_equal(length(unique(res$module_table$module)), res$n_modules)
  expect_false(any(is.na(res$module_table$module)))
  expect_true(all(res$module_table$node %in% igraph::V(g)$name))
  expect_false(any(duplicated(res$module_table$node)))
})

test_that("min_size filters out small modules", {
  set.seed(1)
  loose <- find_modules(g, min_size = 2, plot = FALSE)
  set.seed(1)
  strict <- find_modules(g, min_size = 10, plot = FALSE)

  expect_lte(strict$n_modules, loose$n_modules)

  # Every surviving module is at least min_size nodes.
  sizes <- table(strict$module_table$module)
  if (length(sizes) > 0) expect_gte(min(sizes), 10)
})

test_that("return_subgraphs toggles the subgraphs slot", {
  set.seed(1)
  expect_null(find_modules(g, plot = FALSE, return_subgraphs = FALSE)$subgraphs)

  set.seed(1)
  res <- find_modules(g, plot = FALSE, return_subgraphs = TRUE)
  expect_type(res$subgraphs, "list")
  expect_equal(length(res$subgraphs), res$n_modules)
  for (sg in res$subgraphs) expect_s3_class(sg, "igraph")
})

test_that("several clustering methods are supported", {
  for (m in c("louvain", "walktrap", "fast_greedy", "infomap")) {
    set.seed(1)
    res <- find_modules(g, method = m, plot = FALSE)
    expect_equal(res$method, m)
    expect_gt(res$n_modules, 0)
  }
})

test_that("find_modules() is reproducible under a fixed seed", {
  set.seed(99)
  a <- find_modules(g, method = "louvain", plot = FALSE)
  set.seed(99)
  b <- find_modules(g, method = "louvain", plot = FALSE)

  expect_equal(a$module_table, b$module_table)
  expect_equal(a$n_modules, b$n_modules)
})

test_that("leiden defaults to modularity and does not return one node per community", {
  # igraph's cluster_leiden() defaults to the CPM objective at resolution 1.
  # On a weighted graph whose weights are below 1, no grouping ever beats a
  # singleton under CPM, so every node became its own community and min_size
  # then discarded all of them -- an empty result table with nothing to say why.
  gw <- test_graph_weighted()
  igraph::E(gw)$weight <- igraph::E(gw)$weight / max(igraph::E(gw)$weight)

  res <- find_modules(gw, method = "leiden", plot = FALSE, weights = "weight")
  expect_gt(res$n_modules, 0)
  expect_lt(res$n_modules, igraph::vcount(gw))
  expect_gt(nrow(res$result), 0)

  # The old behaviour is still reachable, and still produces singletons --
  # which is a property of CPM at this resolution, not of netkit.
  cpm <- suppressWarnings(
    find_modules(gw, method = "leiden", objective_function = "CPM",
                 resolution = 1, min_size = 1, plot = FALSE,
                 weights = "weight")
  )
  expect_equal(cpm$n_modules, igraph::vcount(gw))
})

test_that("the leiden controls are honoured and are inert elsewhere", {
  gw <- test_graph_weighted()

  # Higher resolution must not give fewer communities.
  coarse <- find_modules(gw, method = "leiden", resolution = 0.5,
                         min_size = 1, plot = FALSE, weights = "weight")
  fine <- find_modules(gw, method = "leiden", resolution = 2,
                       min_size = 1, plot = FALSE, weights = "weight")
  expect_gte(fine$n_modules, coarse$n_modules)

  expect_error(find_modules(gw, method = "leiden", resolution = 0,
                            plot = FALSE, weights = "weight"),
               "single positive number")

  # Silently ignored arguments make a method look insensitive to tuning.
  expect_warning(find_modules(gw, method = "louvain", resolution = 2,
                              plot = FALSE, weights = "weight"),
                 "method = \"leiden\" only")
})

test_that("find_modules() reports the modularity of the partition it found", {
  g <- test_graph()
  res <- find_modules(g, method = "louvain", plot = FALSE)

  # Must match igraph's own Q for the same membership, computed on the
  # unfiltered partition rather than on whatever survived min_size.
  expect_true(is.finite(res$modularity))
  expect_gt(res$modularity, 0)
  expect_lt(res$modularity, 1)

  # Louvain optimizes Q, so it must beat an arbitrary partition of the graph.
  set.seed(61)
  random_membership <- sample(seq_len(res$n_modules), igraph::vcount(g),
                              replace = TRUE)
  expect_gt(res$modularity, igraph::modularity(g, random_membership))
})

test_that("an empty result says why instead of returning a bare zero-row table", {
  g <- test_graph(n = 40)
  expect_warning(find_modules(g, min_size = igraph::vcount(g) + 1, plot = FALSE),
                 "none has at least")
})
