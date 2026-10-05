g <- test_graph()

test_that("find_modules() returns the documented structure", {
  set.seed(1)
  res <- find_modules(g, plot = FALSE)

  expect_named(res, c("result", "module_table", "n_modules", "subgraphs",
                      "method", "graph"))
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
