# netkit's edge-weight contract. These tests exist because the package used to
# have three different, undocumented weight semantics at once: igraph's implicit
# `weight` attribute was read as a *cost* by the path-based metrics, as a
# *strength* by community detection, and ignored entirely by the diffusion
# matrix code. Attaching a confidence score therefore inverted every path metric
# while leaving diffusion unchanged, and nothing said so.
#
# The central assertion here is the identity property: on a graph whose every
# weight is 1, a weighted call must equal the unweighted one. It holds for every
# metric at once and needs no reference values, so it is the net that catches a
# misaligned weight vector anywhere in the package.

# --- The identity property --------------------------------------------------

test_that("unit weights reproduce the unweighted result exactly", {
  gw <- test_graph_unit_weights()
  gu <- igraph::delete_edge_attr(gw, "weight")

  # summarize_graph_metrics(): Louvain is stochastic, so both runs are seeded.
  # Is_weighted and Avg_strength exist only to make the two row-bindable.
  set.seed(5); a <- summarize_graph_metrics(gw, weights = "weight")
  set.seed(5); b <- summarize_graph_metrics(gu)
  shared <- setdiff(names(a), c("Is_weighted", "Avg_strength"))
  expect_equal(a[shared], b[shared], tolerance = 1e-8)
  expect_true(a$Is_weighted)
  expect_false(b$Is_weighted)

  # find_hubs() / find_bottlenecks()
  for (fn in list(find_hubs, find_bottlenecks)) {
    wa <- fn(gw, plot = FALSE, weights = "weight")$result
    wb <- fn(gu, plot = FALSE)$result
    expect_equal(wa$degree_metric, wb$degree_metric, tolerance = 1e-8)
    expect_equal(wa$betweenness, wb$betweenness, tolerance = 1e-8)
    expect_equal(wa[[ncol(wa)]], wb[[ncol(wb)]])   # the is_* flag
  }

  # network_diffusion(), all three methods
  for (m in c("laplacian", "heat", "rwr")) {
    da <- network_diffusion(gw, c("w1", "w2"), method = m, weights = "weight")
    db <- network_diffusion(gu, c("w1", "w2"), method = m)
    expect_equal(da$score, db$score, tolerance = 1e-10, info = m)
  }

  # calculate_roles(): a fixed membership removes the clustering's randomness,
  # isolating the weighted arithmetic.
  mem <- igraph::membership(igraph::cluster_louvain(gu))
  ra <- calculate_roles(gw, communities = mem, plot = FALSE, weights = "weight")$result
  rb <- calculate_roles(gu, communities = mem, plot = FALSE)$result
  expect_equal(ra$z, rb$z, tolerance = 1e-8)
  expect_equal(ra$p, rb$p, tolerance = 1e-8)
  expect_identical(ra$role, rb$role)

  # robustness_analysis()
  aa <- robustness_analysis(gw, removal_strategy = "degree", steps = 6,
                            plot = FALSE, weights = "weight")
  bb <- robustness_analysis(gu, removal_strategy = "degree", steps = 6,
                            plot = FALSE)
  expect_equal(aa$auc$efficiency, bb$auc$efficiency, tolerance = 1e-8)
  expect_equal(aa$auc$lcc_size, bb$auc$lcc_size, tolerance = 1e-8)
})

test_that("weights that are not all equal do change the answers", {
  # The mirror of the identity test: if nothing responded to real weights, the
  # identity above would pass vacuously -- which is exactly how the silent drop
  # in the diffusion kernel survived.
  gw <- test_graph_weighted()
  gu <- igraph::delete_edge_attr(gw, "weight")

  set.seed(5); a <- summarize_graph_metrics(gw, weights = "weight")
  set.seed(5); b <- summarize_graph_metrics(gu)
  responsive <- c("Diameter", "Average_path_length", "Clustering_coefficient",
                  "Avg_betweenness", "Algebraic_connectivity", "Gini_degree")
  for (m in responsive) {
    expect_false(isTRUE(all.equal(a[[m]], b[[m]])), info = m)
  }

  for (m in c("laplacian", "heat", "rwr")) {
    da <- network_diffusion(gw, c("w1", "w2"), method = m, weights = "weight")
    db <- network_diffusion(gu, c("w1", "w2"), method = m)
    expect_false(isTRUE(all.equal(da$score, db$score)), info = m)
  }
})

# --- Strength versus distance ----------------------------------------------

test_that("the two weight_types are reciprocal, not interchangeable", {
  g <- test_graph_weighted()
  w <- igraph::E(g)$weight

  # Declaring w a strength must equal declaring 1/w a distance: that is the
  # conversion as_netkit_weights() performs, checked from the outside.
  as_strength <- find_hubs(g, plot = FALSE, weights = w,
                           weight_type = "strength")$result
  as_distance <- find_hubs(g, plot = FALSE, weights = 1 / w,
                           weight_type = "distance")$result
  expect_equal(as_strength$betweenness, as_distance$betweenness, tolerance = 1e-8)

  # And the two readings of the *same* numbers must differ, or the declaration
  # would be doing nothing.
  flipped <- find_hubs(g, plot = FALSE, weights = w,
                       weight_type = "distance")$result
  expect_false(isTRUE(all.equal(as_strength$betweenness, flipped$betweenness)))
})

test_that("a strong edge is a short path and a weak edge is a long one", {
  # The sign of the effect, stated as a property rather than a recorded number.
  # Getting this backwards is the defect the contract exists to prevent.
  g <- barbell_graph()
  igraph::E(g)$w <- 1
  g_strong <- igraph::set_edge_attr(g, "w", index = igraph::E(g)["L1" %--% "R1"],
                                    value = 10)
  g_weak <- igraph::set_edge_attr(g, "w", index = igraph::E(g)["L1" %--% "R1"],
                                  value = 0.01)

  d_strong <- summarize_graph_metrics(g_strong, weights = "w")$Average_path_length
  d_weak <- summarize_graph_metrics(g_weak, weights = "w")$Average_path_length

  expect_lt(d_strong, d_weak)
})

# --- Diffusion: the barrier property ---------------------------------------

test_that("a nearly severed bridge blocks diffusion", {
  # This is the test that would have caught weights being dropped from the
  # diffusion kernel: with the bridge at weight 0.001 almost no signal should
  # reach the far clique, whereas the unweighted kernel spreads freely.
  g <- barbell_graph()
  igraph::E(g)$w <- 1
  g <- igraph::set_edge_attr(g, "w", index = igraph::E(g)["L1" %--% "R1"],
                             value = 0.001)
  far <- paste0("R", 1:4)

  open <- network_diffusion(g, "L2", method = "rwr")
  blocked <- network_diffusion(g, "L2", method = "rwr", weights = "w")

  expect_lt(mean(blocked$score[blocked$node %in% far]),
            mean(open$score[open$node %in% far]) / 10)
})

test_that("a zero strength is treated as no connection rather than an error", {
  g <- barbell_graph()
  igraph::E(g)$w <- 1
  g <- igraph::set_edge_attr(g, "w", index = igraph::E(g)["L1" %--% "R1"],
                             value = 0)

  # Zero strength means infinite distance, which the path metrics must handle.
  expect_no_error(summarize_graph_metrics(g, weights = "w"))
  res <- network_diffusion(g, "L2", method = "rwr", weights = "w")
  expect_equal(sum(res$score[res$node %in% paste0("R", 1:4)]), 0,
               tolerance = 1e-12)
})

# --- Kernel reuse ----------------------------------------------------------

test_that("a kernel built with different weights is rejected, not reused", {
  g <- test_graph_weighted()

  weighted_kernel <- prepare_diffusion(g, method = "rwr", weights = "weight")
  # The graph carries a `weight` attribute, so building a deliberately unweighted
  # kernel from it warns -- which is the point of the default, and is asserted in
  # its own test below.
  plain_kernel <- suppressWarnings(prepare_diffusion(g, method = "rwr"))

  expect_error(
    network_diffusion(g, c("w1", "w2"), method = "rwr",
                      precompute = weighted_kernel),
    "built with different edge weights"
  )
  expect_error(
    network_diffusion(g, c("w1", "w2"), method = "rwr", weights = "weight",
                      precompute = plain_kernel),
    "built with different edge weights"
  )
  # Matching weights are accepted.
  expect_no_error(
    network_diffusion(g, c("w1", "w2"), method = "rwr", weights = "weight",
                      precompute = weighted_kernel)
  )
})

# --- The ignore-and-warn default -------------------------------------------

test_that("a weight attribute that is being ignored produces a warning", {
  g <- test_graph_weighted()

  expect_warning(summarize_graph_metrics(g), "being ignored")
  expect_warning(find_hubs(g, plot = FALSE), "being ignored")
  expect_warning(prepare_diffusion(g, method = "rwr"), "being ignored")

  # Naming the attribute, or a graph that has none, is silent.
  expect_no_warning(summarize_graph_metrics(g, weights = "weight"))
  expect_no_warning(summarize_graph_metrics(igraph::delete_edge_attr(g, "weight")))
})

test_that("an edge attribute under any other name is not picked up implicitly", {
  # Only `weight` triggers the warning, because only `weight` is what igraph
  # would have read behind netkit's back.
  g <- test_graph_weighted()
  g <- igraph::set_edge_attr(igraph::delete_edge_attr(g, "weight"),
                             "confidence", value = runif(igraph::ecount(g), .1, 1))
  expect_no_warning(summarize_graph_metrics(g))
  expect_no_warning(summarize_graph_metrics(g, weights = "confidence"))
})

# --- Validation ------------------------------------------------------------

test_that("malformed weights are rejected eagerly with a specific message", {
  g <- test_graph_gnp()
  ne <- igraph::ecount(g)

  expect_error(summarize_graph_metrics(g, weights = "nope"),
               "Edge attribute 'nope' not found")
  expect_error(summarize_graph_metrics(g, weights = rep(1, ne - 1)),
               "has length")
  expect_error(summarize_graph_metrics(g, weights = c(rep(1, ne - 1), NA)),
               "missing values")
  expect_error(summarize_graph_metrics(g, weights = c(rep(1, ne - 1), Inf)),
               "non-finite")
  expect_error(summarize_graph_metrics(g, weights = c(rep(1, ne - 1), -1)),
               "negative values")
  expect_error(summarize_graph_metrics(g, weights = rep(0, ne)),
               "zero for every edge")
  expect_error(summarize_graph_metrics(g, weights = TRUE),
               "must be NULL, the name of an edge attribute")
})

test_that("a signed attribute is refused rather than silently misused", {
  # Activation/inhibition and signed correlations are not edge weights. Taking
  # them as such would make the reciprocal conversion meaningless.
  g <- test_graph_gnp()
  signed <- rep(c(-1, 1), length.out = igraph::ecount(g))
  expect_error(summarize_graph_metrics(g, weights = signed), "negative values")
})

# --- Weights survive the transformations that renumber edges ---------------

test_that("collapsing a directed graph keeps weights aligned with edges", {
  # as_undirected(mode = "collapse") merges reciprocal pairs and renumbers the
  # edges, so a positional weight vector would silently refer to the wrong ones.
  # Reciprocal strengths are summed, so a graph with every reciprocal pair
  # present has twice the total strength of its undirected form.
  g <- igraph::make_full_graph(5, directed = TRUE)
  igraph::V(g)$name <- paste0("v", 1:5)
  igraph::E(g)$w <- 1

  m <- summarize_graph_metrics(g, weights = "w")
  expect_true(m$Is_weighted)
  expect_false(m$Is_directed == FALSE)   # reported as having been directed
  # Every pair is reciprocal, so each collapsed edge has strength 2.
  expect_equal(m$Avg_strength, 2 * (igraph::vcount(g) - 1))
})

test_that("spinglass on a disconnected graph keeps its subgraph weights aligned", {
  # find_modules() restricts spinglass to the largest component, which drops
  # edges. The weights must be subset with them.
  g <- test_graph_disconnected()
  igraph::E(g)$w <- seq_len(igraph::ecount(g))

  set.seed(1)
  expect_warning(
    find_modules(g, method = "spinglass", min_size = 1, plot = FALSE,
                 weights = "w"),
    "cannot work with unconnected graph"
  )

  set.seed(1)
  res <- suppressWarnings(
    find_modules(g, method = "spinglass", min_size = 1, plot = FALSE,
                 weights = "w")
  )
  expect_s3_class(res$result, "data.frame")
  expect_gt(res$n_modules, 0)
  # Only the largest component is covered, which is the documented behaviour.
  expect_lt(nrow(res$result), igraph::vcount(g))
})

test_that("fluid_communities warns that it cannot use weights", {
  g <- test_graph_gnp()
  igraph::E(g)$w <- runif(igraph::ecount(g), 0.1, 1)
  set.seed(1)
  expect_warning(
    find_modules(g, method = "fluid_communities", no.of.communities = 3,
                 plot = FALSE, weights = "w"),
    "does not support edge weights"
  )
})

# --- edge_betweenness is the algorithm that needs costs --------------------

test_that("edge_betweenness clustering receives costs, not strengths", {
  # Every other community algorithm reads a weight as a strength; this one
  # removes high-betweenness edges and so needs the cost vector. If it were
  # handed strengths, declaring the same numbers a distance instead would be a
  # no-op -- this asserts the two readings actually differ.
  g <- test_graph_gnp(n = 30, p = 0.12, seed = 21)
  igraph::E(g)$w <- runif(igraph::ecount(g), 0.1, 1)

  # igraph itself notes that a weighted run selects membership by modularity;
  # asserted rather than ignored, because an unasserted warning is a signal.
  set.seed(1)
  expect_warning(
    find_modules(g, method = "edge_betweenness", min_size = 1, plot = FALSE,
                 weights = "w", weight_type = "strength"),
    "highest modularity score"
  )

  set.seed(1)
  as_s <- suppressWarnings(find_modules(g, method = "edge_betweenness", min_size = 1,
                       plot = FALSE, weights = "w",
                       weight_type = "strength")$result)
  set.seed(1)
  as_d <- suppressWarnings(find_modules(g, method = "edge_betweenness", min_size = 1,
                       plot = FALSE, weights = "w",
                       weight_type = "distance")$result)

  expect_false(isTRUE(all.equal(as_s$module, as_d$module)))
})

# --- Seed weights ----------------------------------------------------------

test_that("seed_weights scales the diffusion linearly", {
  g <- test_graph_gnp()
  seeds <- c("m1", "m2", "m3")

  base <- network_diffusion(g, seeds, method = "rwr")
  doubled <- network_diffusion(g, seeds, method = "rwr",
                               seed_weights = c(2, 2, 2))

  # All three kernels are linear operators on f0, so scaling every seed by 2
  # must scale every score by 2. The tolerance is the RWR convergence threshold,
  # which is relative to sum(abs(f0)) precisely so that this holds at any scale.
  expect_equal(doubled$score[match(base$node, doubled$node)], 2 * base$score,
               tolerance = 1e-5)
})

test_that("seed_weights can be named or positional and are validated", {
  g <- test_graph_gnp()
  seeds <- c("m1", "m2", "m3")

  named <- network_diffusion(g, seeds, method = "rwr",
                             seed_weights = c(m1 = 3, m2 = 1, m3 = 2))
  positional <- network_diffusion(g, seeds, method = "rwr",
                                  seed_weights = c(3, 1, 2))
  expect_equal(named$score, positional$score, tolerance = 1e-10)

  expect_error(network_diffusion(g, seeds, seed_weights = c(1, 2)),
               "length 2 but there are 3")
  expect_error(network_diffusion(g, seeds, seed_weights = c(m1 = 1, m2 = 2)),
               "no entry for")
  expect_error(network_diffusion(g, seeds, seed_weights = c(1, 2, NA)),
               "missing or non-finite")
})

test_that("unknown seed nodes warn, and all-unknown errors", {
  g <- test_graph_gnp()
  expect_warning(network_diffusion(g, c("m1", "ghost"), method = "rwr"),
                 "not vertices of the graph")
  expect_error(network_diffusion(g, c("ghost1", "ghost2"), method = "rwr"),
               "None of the 'seed_nodes'")
  expect_error(network_diffusion(g, character(0), method = "rwr"),
               "needs at least one seed")
})

# --- robustness_analysis's new strategy ------------------------------------

test_that("removal_strategy = 'strength' requires weights", {
  g <- test_graph_weighted()
  expect_error(
    robustness_analysis(igraph::delete_edge_attr(g, "weight"),
                        removal_strategy = "strength", steps = 5, plot = FALSE),
    "needs edge weights"
  )
  expect_no_error(
    robustness_analysis(g, removal_strategy = "strength", steps = 5,
                        plot = FALSE, weights = "weight")
  )
})

# --- plot_CCDF's strength path ---------------------------------------------

test_that("plot_CCDF plots a strength distribution when weighted", {
  g <- test_graph_weighted()

  p_deg <- plot_CCDF(igraph::delete_edge_attr(g, "weight"), show_PL = FALSE)
  p_str <- plot_CCDF(g, show_PL = FALSE, weights = "weight")

  expect_s3_class(p_str, "ggplot")
  expect_equal(p_str$labels$x, "Strength, s")
  expect_equal(p_deg$labels$x, "Degree, k")

  # A strength CCDF is evaluated on continuous values, so it is not the integer
  # degree sequence.
  expect_false(all(p_str$data$degree == round(p_str$data$degree)))

  # It is still a proper CCDF: non-increasing and bounded by 1.
  expect_true(all(diff(p_str$data$ccdf) <= 1e-12))
  expect_lte(max(p_str$data$ccdf), 1)
})
