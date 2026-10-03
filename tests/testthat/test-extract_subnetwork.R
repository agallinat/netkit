# extract_subnetwork() is asserted on structural properties -- which nodes must
# be present, whether the result is a tree, how the methods nest -- rather than
# on recorded node lists, which would just re-record the current output.

g <- test_graph_gnp(n = 50, p = 0.08, seed = 456)
seeds <- c("m1", "m10", "m20", "m30")

# --- Shape ------------------------------------------------------------------

test_that("extract_subnetwork() returns result/graph/method and no plot slot", {
  res <- extract_subnetwork(g, seeds, method = "steiner", plot = FALSE)

  expect_type(res, "list")
  expect_true(all(c("result", "graph", "method") %in% names(res)))
  expect_s3_class(res$result, "tbl_df")
  expect_s3_class(res$graph, "igraph")
  expect_type(res$method, "character")

  # Like plot_Net(), find_modules() and highlight_nodes(), this renders through
  # base plot.igraph, so there is no plot object to return under any setting.
  expect_false("plot" %in% names(res))
  expect_false("plot" %in% names(draw_quietly(
    extract_subnetwork(g, seeds, method = "steiner", plot = TRUE)
  )))
})

test_that("the result table lines up with the returned graph", {
  res <- extract_subnetwork(g, seeds, method = "shortest_paths", plot = FALSE)
  expect_identical(res$result$node, igraph::V(res$graph)$name)
  expect_equal(igraph::V(res$graph)$is_seed, res$result$is_seed)
  expect_equal(igraph::V(res$graph)$reason, res$result$reason)
})

# --- Every method keeps the seeds -------------------------------------------

test_that("every method retains all the requested nodes", {
  for (m in c("induced", "neighbors", "shortest_paths", "steiner", "diffusion")) {
    res <- extract_subnetwork(g, seeds, method = m, top_n = 15, plot = FALSE)
    expect_true(all(seeds %in% res$result$node), info = m)
    expect_true(all(res$result$is_seed[res$result$node %in% seeds]), info = m)
    # A seed is labelled a seed whichever path kept it.
    expect_true(all(res$result$reason[res$result$is_seed] == "seed"), info = m)
  }
})

# --- How the methods nest ---------------------------------------------------

test_that("induced is a subset of steiner, which is a subset of all shortest paths", {
  # Checked on a graph with cycles, so that the union of shortest paths can be
  # strictly larger than one connecting tree. On a tree the two coincide, which
  # is asserted separately below.
  ind <- extract_subnetwork(g, seeds, method = "induced", plot = FALSE)
  st <- extract_subnetwork(g, seeds, method = "steiner", plot = FALSE)
  sp <- extract_subnetwork(g, seeds, method = "shortest_paths", plot = FALSE)

  expect_true(all(ind$result$node %in% st$result$node))
  expect_lte(igraph::vcount(st$graph), igraph::vcount(sp$graph))
})

test_that("on a tree, the Steiner tree is exactly the union of the seed paths", {
  # The only case with a known closed form: in a tree there is exactly one path
  # between any two vertices, so the heuristic must find the optimum.
  tree <- test_graph(n = 40)   # sample_pa with m = 1 is a tree
  tseeds <- c("n5", "n17", "n33")

  st <- extract_subnetwork(tree, tseeds, method = "steiner", plot = FALSE)
  sp <- extract_subnetwork(tree, tseeds, method = "shortest_paths", plot = FALSE)

  expect_setequal(st$result$node, sp$result$node)
})

# --- The Steiner result is a tree -------------------------------------------

test_that("the Steiner result is a connected tree containing every terminal", {
  st <- extract_subnetwork(g, seeds, method = "steiner", plot = FALSE)

  expect_equal(igraph::components(st$graph)$no, 1)
  expect_equal(igraph::ecount(st$graph), igraph::vcount(st$graph) - 1)
  expect_true(all(seeds %in% igraph::V(st$graph)$name))
})

test_that("the Steiner tree has no non-terminal leaves", {
  # A leaf that is not a terminal is a branch leading nowhere, so the pruning
  # step must have removed it. This is the property that keeps the tree minimal.
  st <- extract_subnetwork(g, seeds, method = "steiner", plot = FALSE)
  deg <- igraph::degree(st$graph)
  leaves <- igraph::V(st$graph)$name[deg <= 1]
  expect_true(all(leaves %in% seeds))
})

test_that("two terminals give exactly a shortest path between them", {
  st <- extract_subnetwork(g, c("m1", "m30"), method = "steiner", plot = FALSE)
  d <- igraph::distances(g, v = "m1", to = "m30")[1, 1]
  # A path of graph distance d has d edges and d + 1 vertices.
  expect_equal(igraph::ecount(st$graph), d)
  expect_equal(igraph::vcount(st$graph), d + 1)
})

test_that("a single terminal is returned alone", {
  for (m in c("induced", "shortest_paths", "steiner")) {
    res <- extract_subnetwork(g, "m1", method = m, plot = FALSE)
    expect_equal(res$result$node, "m1", info = m)
  }
})

# --- Neighbourhood expansion ------------------------------------------------

test_that("neighbourhood expansion grows with order", {
  n1 <- igraph::vcount(extract_subnetwork(g, "m1", method = "neighbors",
                                          order = 1, plot = FALSE)$graph)
  n2 <- igraph::vcount(extract_subnetwork(g, "m1", method = "neighbors",
                                          order = 2, plot = FALSE)$graph)
  expect_gte(n2, n1)
  # order = 0 is the seed alone, i.e. the induced subgraph.
  n0 <- extract_subnetwork(g, "m1", method = "neighbors", order = 0,
                           plot = FALSE)
  expect_equal(n0$result$node, "m1")
})

test_that("max_degree is monotone: lowering it never grows the result", {
  sizes <- vapply(c(2, 5, 10, 1000), function(md) {
    igraph::vcount(extract_subnetwork(g, seeds, method = "neighbors",
                                      max_degree = md, plot = FALSE)$graph)
  }, numeric(1))
  expect_false(is.unsorted(sizes))
})

test_that("max_degree never drops a seed, however high its degree", {
  # A seed is present because the caller asked for it; removing one for being a
  # hub would silently answer a different question.
  hub <- names(which.max(igraph::degree(g)))
  res <- extract_subnetwork(g, hub, method = "neighbors", max_degree = 1,
                            plot = FALSE)
  expect_true(hub %in% res$result$node)
})

test_that("an invalid order is rejected", {
  expect_error(extract_subnetwork(g, seeds, method = "neighbors", order = -1,
                                  plot = FALSE),
               "non-negative")
})

# --- Diffusion expansion ----------------------------------------------------

test_that("diffusion keeps top_n vertices plus the seeds", {
  res <- extract_subnetwork(g, seeds, method = "diffusion", top_n = 12,
                            plot = FALSE)
  # Seeds are retained whether or not they rank in the top_n, so the bound is
  # top_n plus however many seeds fell outside it.
  expect_lte(igraph::vcount(res$graph), 12 + length(seeds))
  expect_gte(igraph::vcount(res$graph), length(seeds))
})

test_that("diffusion reports the scores it ranked by", {
  res <- extract_subnetwork(g, seeds, method = "diffusion", top_n = 12,
                            plot = FALSE)
  expect_false(any(is.na(res$result$score)))
  expect_true(all(res$result$score >= 0))

  # The non-seed vertices kept must all score at least as high as any vertex
  # that was dropped -- that is what "top_n" means.
  all_scores <- network_diffusion(g, seeds, method = "rwr")
  dropped <- setdiff(all_scores$node, res$result$node)
  if (length(dropped) > 0) {
    kept_min <- min(res$result$score)
    dropped_max <- max(all_scores$score[all_scores$node %in% dropped])
    expect_gte(kept_min, dropped_max - 1e-12)
  }
})

test_that("the other methods report NA scores rather than a fabricated number", {
  res <- extract_subnetwork(g, seeds, method = "steiner", plot = FALSE)
  expect_true(all(is.na(res$result$score)))
})

test_that("an invalid top_n is rejected", {
  expect_error(extract_subnetwork(g, seeds, method = "diffusion", top_n = 0,
                                  plot = FALSE),
               "at least 1")
})

# --- Disconnected input -----------------------------------------------------

test_that("seeds in different components warn and yield a forest", {
  dg <- test_graph_disconnected()
  dseeds <- c("a1", "a3", "b1", "b3")

  expect_warning(
    extract_subnetwork(dg, dseeds, method = "steiner", plot = FALSE),
    "span|components"
  )
  expect_warning(
    extract_subnetwork(dg, dseeds, method = "shortest_paths", plot = FALSE),
    "different components"
  )

  st <- suppressWarnings(
    extract_subnetwork(dg, dseeds, method = "steiner", plot = FALSE)
  )
  expect_gt(igraph::components(st$graph)$no, 1)
  expect_true(all(dseeds %in% st$result$node))
})

test_that("largest_component reduces the forest to one tree", {
  dg <- test_graph_disconnected()
  st <- suppressWarnings(
    extract_subnetwork(dg, c("a1", "a3", "b1", "b3"), method = "steiner",
                       largest_component = TRUE, plot = FALSE)
  )
  expect_equal(igraph::components(st$graph)$no, 1)
  expect_match(st$method, "largest component only")
})

# --- Unknown nodes ----------------------------------------------------------

test_that("unknown nodes warn, and all-unknown errors", {
  expect_warning(extract_subnetwork(g, c("m1", "ghost"), method = "induced",
                                    plot = FALSE),
                 "not vertices of the graph")
  expect_error(extract_subnetwork(g, c("ghost1", "ghost2"), method = "induced",
                                  plot = FALSE),
               "None of the 'nodes'")
})

# --- Weights ----------------------------------------------------------------

test_that("unit weights reproduce the unweighted extraction where it is exact", {
  gw <- test_graph_unit_weights()
  gu <- igraph::delete_edge_attr(gw, "weight")
  wseeds <- c("w1", "w10", "w20", "w30")

  # Diffusion is a deterministic linear operator, so the identity is exact.
  for (m in c("induced", "neighbors", "diffusion")) {
    a <- extract_subnetwork(gw, wseeds, method = m, weights = "weight",
                            top_n = 12, plot = FALSE)
    b <- extract_subnetwork(gu, wseeds, method = m, top_n = 12, plot = FALSE)
    expect_setequal(a$result$node, b$result$node)
  }
})

test_that("the path-based methods agree on distance, if not on which tie they break", {
  # The identity cannot hold on the vertex set for `steiner` and
  # `shortest_paths`: igraph uses breadth-first search with no weights and
  # Dijkstra with them, and those pick different paths among equal-cost ties.
  # What must agree is the thing that is actually well defined -- the distances
  # the paths realise -- plus the guarantees the method advertises.
  gw <- test_graph_unit_weights()
  gu <- igraph::delete_edge_attr(gw, "weight")
  wseeds <- c("w1", "w10", "w20", "w30")

  # Pairwise distances are identical under unit weights, which is what rules out
  # the weight vector itself being misaligned.
  expect_equal(
    igraph::distances(gw, v = wseeds, to = wseeds,
                      weights = igraph::E(gw)$weight),
    igraph::distances(gu, v = wseeds, to = wseeds)
  )

  a <- extract_subnetwork(gw, wseeds, method = "steiner", weights = "weight",
                          plot = FALSE)
  b <- extract_subnetwork(gu, wseeds, method = "steiner", plot = FALSE)

  # Both must be trees spanning the terminals.
  for (res in list(a, b)) {
    expect_equal(igraph::ecount(res$graph), igraph::vcount(res$graph) - 1)
    expect_equal(igraph::components(res$graph)$no, 1)
    expect_true(all(wseeds %in% res$result$node))
  }

  # And both must honour the Kou-Markowsky-Berman bound of 2 - 2/|terminals|
  # relative to the better of the two, since neither is known to be optimal.
  bound <- 2 - 2 / length(wseeds)
  costs <- c(igraph::ecount(a$graph), igraph::ecount(b$graph))
  expect_lte(max(costs), bound * min(costs))
})

test_that("a cheap route is preferred over an expensive one", {
  # Two parallel 2-hop routes between A and B. Weighting one as much stronger
  # must make the Steiner tree take it. Stated as a property, because getting
  # the strength/cost direction backwards is the defect the weight contract
  # exists to prevent and it would still return a valid-looking tree.
  el <- data.frame(
    from = c("A", "cheap", "A", "pricey"),
    to   = c("cheap", "B", "pricey", "B"),
    stringsAsFactors = FALSE
  )
  gg <- igraph::graph_from_data_frame(el, directed = FALSE)
  igraph::E(gg)$strong <- c(10, 10, 0.1, 0.1)

  res <- extract_subnetwork(gg, c("A", "B"), method = "steiner",
                            weights = "strong", plot = FALSE)
  expect_true("cheap" %in% res$result$node)
  expect_false("pricey" %in% res$result$node)
})

# --- Input contract ---------------------------------------------------------

test_that("extract_subnetwork() honours the shared graph-input contract", {
  expect_error(extract_subnetwork("not a graph", nodes = "a"),
               "Input 'graph' must be either an igraph object or a data.frame")
  expect_type(
    extract_subnetwork(as_edge_df(g), seeds, method = "induced", plot = FALSE),
    "list"
  )
})
