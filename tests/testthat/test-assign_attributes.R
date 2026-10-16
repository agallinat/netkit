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

  # 2 of the ring's 5 edges are covered, which the function warns about. The
  # warning used to be unreachable: it compared nrow(merged), which all.x = TRUE
  # pins to the graph's edge count, against the table's row count.
  expect_warning(g2 <- assign_attributes(g, edge_table = edges),
                 "2 of 5 graph edges matched")
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

  expect_warning(g2 <- assign_attributes(g, nodes_table = nodes, edge_table = edges),
                 "1 of 5 graph edges matched")
  expect_true("grp" %in% igraph::vertex_attr_names(g2))
  expect_true("w" %in% igraph::edge_attr_names(g2))
})

test_that("parallel edges each take their own edge_table row", {
  # An endpoint pair does not identify an edge in a multigraph. merge() formed
  # the cartesian product of the two duplicate groups, and `values[index] <- `
  # then kept whichever row the merge put last, so both parallel edges got one
  # value and which one depended on the table's row order.
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("a", "a", "b"), to = c("b", "b", "c")), directed = TRUE)
  edges <- data.frame(from = c("a", "a", "b"), to = c("b", "b", "c"),
                      label = c("first", "second", "third"))

  # A complete 1:1 pairing, so nothing is unmatched and nothing warns.
  g2 <- assign_attributes(g, edge_table = edges)
  expect_equal(igraph::E(g2)$label, c("first", "second", "third"))

  # The k-th edge of a pair takes the k-th row of that pair, so permuting the
  # table's rows permutes the values rather than collapsing them onto one.
  g3 <- assign_attributes(g, edge_table = edges[c(2, 1, 3), ])
  expect_equal(igraph::E(g3)$label, c("second", "first", "third"))
  # The defect's signature: every parallel edge sharing a single value.
  expect_false(igraph::E(g3)$label[[1]] == igraph::E(g3)$label[[2]])
})

test_that("a pair repeated more often in the graph than in the table warns", {
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("a", "a"), to = c("b", "b")), directed = TRUE)
  edges <- data.frame(from = "a", to = "b", w = 7)

  expect_warning(g2 <- assign_attributes(g, edge_table = edges),
                 "1 of 2 graph edges matched")
  # One edge is covered, the leftover is NA rather than a recycled 7.
  expect_equal(igraph::E(g2)$w, c(7, NA))
})

test_that("an undirected edge_table may list either endpoint first", {
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("a", "b"), to = c("b", "c")), directed = FALSE)
  # Both rows reversed relative to the graph's edge list.
  edges <- data.frame(from = c("b", "c"), to = c("a", "b"), w = c(10, 20))

  g2 <- assign_attributes(g, edge_table = edges)
  expect_equal(igraph::E(g2)$w, c(10, 20))
})

test_that("the unmatched counts share a denominator with what they count", {
  # Reporting a count of matched edges against the table's row count let
  # assign_attributes() claim "5 rows matched" for a 3-row table. Both counts now
  # come from one matching, so each is bounded by its own total.
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("a", "b"), to = c("b", "c")), directed = TRUE)
  edges <- data.frame(from = c("a", "x", "y"), to = c("b", "z", "w"),
                      w = c(1, 2, 3))

  warns <- character()
  withCallingHandlers(
    assign_attributes(g, edge_table = edges),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  # 1 of 3 rows used, 1 of 2 edges matched -- never "3 of 2" or "5 of 3".
  expect_true(any(grepl("1 of 3 rows in 'edge_table'", warns)))
  expect_true(any(grepl("1 of 2 graph edges", warns)))
})

test_that("a node table repeating an id is not reported as unmatched", {
  # match() resolves every repeat to the first row, so counting matched vertices
  # and reporting it against nrow(nodes_table) overstated the unmatched rows.
  g <- ring_graph()
  nodes <- data.frame(node = rep(letters[1:5], each = 2),
                      grp = rep(letters[1:5], each = 2))

  # All 5 vertices match and all 5 distinct ids are present. The repeats are
  # genuinely unused rows, so 5 of 10 is the honest count.
  expect_warning(g2 <- assign_attributes(g, nodes_table = nodes),
                 "5 of 10 rows in 'nodes_table'")
  expect_equal(igraph::V(g2)$grp, letters[1:5])
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
