test_that("node attributes are attached by name", {
  g <- ring_graph()
  nodes <- data.frame(node = c("a", "c", "e"), score2 = c(1.5, 2.3, 3.1))

  # Only 3 of 5 vertices match, which the function warns about.
  expect_warning(g2 <- assign_attributes(g, nodes_table = nodes),
                 "3 of 5 graph vertices matched")

  expect_s3_class(g2, "igraph")
  expect_true("score2" %in% igraph::vertex_attr_names(g2))

  got <- stats::setNames(igraph::V(g2)$score2, igraph::V(g2)$name)
  expect_equal(got[["a"]], 1.5)
  expect_equal(got[["c"]], 2.3)
  expect_equal(got[["e"]], 3.1)
  # Unmatched vertices are left missing rather than silently zeroed.
  expect_true(is.na(got[["b"]]))
  expect_true(is.na(got[["d"]]))
})

test_that("a fully matching node table produces no warning", {
  g <- ring_graph()
  nodes <- data.frame(node = letters[1:5], grp = letters[1:5])

  g2 <- assign_attributes(g, nodes_table = nodes)
  expect_equal(igraph::V(g2)$grp, letters[1:5])
})

test_that("edge attributes are attached by endpoint pair", {
  g <- ring_graph()
  edges <- data.frame(from = c("a", "b"), to = c("b", "c"), weight2 = c(10, 20))

  g2 <- assign_attributes(g, edge_table = edges)
  expect_true("weight2" %in% igraph::edge_attr_names(g2))

  ends <- igraph::ends(g2, igraph::E(g2), names = TRUE)
  w <- igraph::E(g2)$weight2
  ab <- which((ends[, 1] == "a" & ends[, 2] == "b") |
                (ends[, 1] == "b" & ends[, 2] == "a"))
  expect_equal(w[ab], 10)
})

test_that("node and edge tables can be supplied together", {
  g <- ring_graph()
  nodes <- data.frame(node = letters[1:5], grp = 1:5)
  edges <- data.frame(from = "a", to = "b", w = 99)

  g2 <- assign_attributes(g, nodes_table = nodes, edge_table = edges)
  expect_true("grp" %in% igraph::vertex_attr_names(g2))
  expect_true("w" %in% igraph::edge_attr_names(g2))
})

test_that("the graph is returned unchanged when no tables are supplied", {
  g <- ring_graph()
  expect_equal(igraph::vertex_attr_names(assign_attributes(g)),
               igraph::vertex_attr_names(g))
})

test_that("overwrite = FALSE preserves an existing attribute", {
  g <- ring_graph()
  igraph::V(g)$grp <- rep("original", 5)
  nodes <- data.frame(node = letters[1:5], grp = rep("new", 5))

  kept <- assign_attributes(g, nodes_table = nodes, overwrite = FALSE)
  expect_equal(igraph::V(kept)$grp, rep("original", 5))

  # overwrite = TRUE announces itself; assert that rather than letting the warning
  # sit in the suite's summary as ambient noise.
  expect_warning(replaced <- assign_attributes(g, nodes_table = nodes, overwrite = TRUE),
                 "Overwriting existing vertex attribute 'grp'")
  expect_equal(igraph::V(replaced)$grp, rep("new", 5))
})

test_that("a diffusion result can be attached directly to the graph", {
  # The documented chaining idiom: diffusion scores become a vertex attribute.
  g <- test_graph()
  diffusion <- network_diffusion(g, seed_nodes = c("n1", "n2"), method = "rwr")
  names(diffusion)[1] <- "name"

  # The fixture already carries a numeric `score`, so this overwrites it.
  expect_warning(g2 <- assign_attributes(g, nodes_table = diffusion),
                 "Overwriting existing vertex attribute 'score'")
  expect_true("score" %in% igraph::vertex_attr_names(g2))
  expect_false(any(is.na(igraph::V(g2)$score)))
})
