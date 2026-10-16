# Null models and the significance tests built on them. The sharpest available
# invariant is that a degree-preserving model preserves the degree sequence
# *exactly*, so most of these assert that rather than any distributional claim.

g <- test_graph()

# --- null_model() -----------------------------------------------------------

test_that("null_model() returns n graphs with the recorded model", {
  for (m in c("rewire", "configuration", "erdos_renyi")) {
    nulls <- null_model(g, model = m, n = 5, seed = 1)
    expect_length(nulls, 5)
    expect_s3_class(nulls, "netkit_null")
    expect_equal(attr(nulls, "model"), m)
    expect_true(all(vapply(nulls, igraph::is_igraph, logical(1))))
  }
})

test_that("the degree-preserving models preserve the degree sequence exactly", {
  # Not "approximately" and not "on average": exactly. This is the one assertion
  # that cannot pass by accident if the model is wrong.
  observed <- sort(unname(igraph::degree(g)))
  for (m in c("rewire", "configuration")) {
    nulls <- null_model(g, model = m, n = 10, seed = 1)
    for (i in seq_along(nulls)) {
      expect_identical(sort(unname(igraph::degree(nulls[[i]]))), observed,
                       info = paste(m, i))
    }
  }
})

test_that("erdos_renyi preserves only the vertex and edge counts", {
  nulls <- null_model(g, model = "erdos_renyi", n = 5, seed = 1)
  for (null_g in nulls) {
    expect_equal(igraph::vcount(null_g), igraph::vcount(g))
    expect_equal(igraph::ecount(null_g), igraph::ecount(g))
  }
  # And it must actually destroy the degree sequence, or it would be the same
  # null as the others under a different name.
  expect_false(identical(sort(unname(igraph::degree(nulls[[1]]))),
                         sort(unname(igraph::degree(g)))))
})

test_that("the ensemble is not n copies of one graph", {
  # A rewiring that never fired would preserve the degree sequence perfectly and
  # pass the test above while being useless as a null.
  nulls <- null_model(g, model = "rewire", n = 5, seed = 1)
  edge_sets <- vapply(nulls, function(x) {
    paste(sort(apply(igraph::as_edgelist(x), 1, paste, collapse = "-")),
          collapse = "|")
  }, character(1))
  expect_gt(length(unique(edge_sets)), 1)
  # And the nulls must differ from the observed graph too.
  observed_edges <- paste(sort(apply(igraph::as_edgelist(g), 1, paste,
                                     collapse = "-")), collapse = "|")
  expect_false(all(edge_sets == observed_edges))
})

test_that("a seed makes the ensemble reproducible and does not leak", {
  a <- null_model(g, model = "rewire", n = 3, seed = 7)
  b <- null_model(g, model = "rewire", n = 3, seed = 7)
  expect_equal(igraph::as_edgelist(a[[1]]), igraph::as_edgelist(b[[1]]))

  # seed = NULL must not call set.seed(NULL), which would re-seed from the clock
  # and destroy the caller's stream -- the defect robustness_analysis() had.
  set.seed(42); invisible(null_model(g, n = 2)); first <- runif(3)
  set.seed(42); invisible(null_model(g, n = 2)); second <- runif(3)
  expect_equal(first, second)
})

test_that("vertex names are carried onto the null graphs", {
  # Needed by anything that keys results by vertex name, which is most of netkit.
  for (m in c("rewire", "configuration", "erdos_renyi")) {
    nulls <- null_model(g, model = m, n = 2, seed = 1)
    expect_setequal(igraph::V(nulls[[1]])$name, igraph::V(g)$name)
  }
})

test_that("weights are permuted, preserving the distribution but not the placement", {
  gw <- test_graph_weighted()
  nulls <- null_model(gw, model = "rewire", n = 5, seed = 1)

  for (null_g in nulls) {
    expect_true("weight" %in% igraph::edge_attr_names(null_g))
    # Same multiset of weights, so any weighted metric is tested against where
    # the weights sit rather than against what they are.
    expect_equal(sort(igraph::E(null_g)$weight), sort(igraph::E(gw)$weight))
  }

  expect_false("weight" %in% igraph::edge_attr_names(
    null_model(gw, n = 1, shuffle_weights = FALSE, seed = 1)[[1]]
  ))
})

test_that("degenerate and invalid input is rejected", {
  empty <- igraph::make_empty_graph(n = 5, directed = FALSE)
  expect_error(null_model(empty), "no edges")
  expect_error(null_model(g, n = 0), "at least 1")
  expect_error(null_model("not a graph"),
               "Input 'graph' must be either an igraph object or a data.frame")
})

test_that("a directed graph is collapsed with a message", {
  expect_message(null_model(test_graph_directed(), n = 2, seed = 1),
                 "converted to undirected")
})

test_that("configuration falls back when vl cannot realize the degree sequence", {
  # "vl" needs a connected realization to exist, which isolated vertices rule
  # out. The fallback must warn rather than silently change model.
  iso <- igraph::add_vertices(test_graph(n = 30), 3,
                              name = c("iso1", "iso2", "iso3"))
  expect_warning(null_model(iso, model = "configuration", n = 2, seed = 1),
                 "configuration.simple")
})

# --- metric_significance() --------------------------------------------------

test_that("metric_significance() returns the full plot/result/graph/method shape", {
  res <- metric_significance(g, metrics = "Modularity", n = 10, seed = 1,
                             plot = FALSE)
  expect_true(all(c("plot", "result", "graph", "method") %in% names(res)))
  expect_null(res$plot)
  expect_s3_class(res$result, "tbl_df")
  expect_s3_class(res$graph, "igraph")
  expect_named(res$result, c("metric", "observed", "null_mean", "null_sd", "z",
                             "p_empirical", "ci_lower", "ci_upper"))
})

test_that("empirical p-values respect the (r+1)/(n+1) floor", {
  n_null <- 20
  res <- metric_significance(g, metrics = c("Modularity",
                                            "Clustering_coefficient"),
                             n = n_null, seed = 1, plot = FALSE)
  # Never zero, and never below the smallest attainable value. Reporting p = 0
  # from a permutation test is a category error: it claims more evidence than
  # the number of permutations can supply.
  expect_true(all(res$result$p_empirical >= 1 / (n_null + 1) - 1e-12))
  expect_true(all(res$result$p_empirical <= 1))
  expect_false(any(res$result$p_empirical == 0))
  expect_match(res$method, "smallest attainable p-value")
})

test_that("the observed value lies inside the reported interval ordering", {
  res <- metric_significance(g, metrics = c("Modularity",
                                            "Clustering_coefficient"),
                             n = 20, seed = 1, plot = FALSE)
  expect_true(all(res$result$ci_lower <= res$result$ci_upper))
  expect_true(all(res$result$null_sd >= 0))
})

test_that("a graph drawn from the null is not significant against it", {
  # The calibration check: an Erdos-Renyi graph tested against an Erdos-Renyi
  # null should show no effect. If this failed, every p-value the function
  # produces would be suspect.
  set.seed(11)
  er <- igraph::sample_gnp(80, 0.08, directed = FALSE)
  res <- metric_significance(er, metrics = "Clustering_coefficient",
                             model = "erdos_renyi", n = 50, seed = 2,
                             plot = FALSE)
  expect_lt(abs(res$result$z), 3)
  expect_gt(res$result$p_empirical, 0.05)
})

test_that("a modular graph is significantly modular against a random null", {
  # The power check, the mirror of the calibration one: three near-cliques with
  # sparse interconnection must beat an Erdos-Renyi null on modularity.
  set.seed(12)
  blocks <- igraph::sample_islands(islands.n = 3, islands.size = 20,
                                   islands.pin = 0.5, n.inter = 2)
  res <- metric_significance(blocks, metrics = "Modularity",
                             model = "erdos_renyi", n = 30, seed = 3,
                             plot = FALSE)
  expect_gt(res$result$z, 3)
  expect_lte(res$result$p_empirical, 1 / 31 + 1e-12)
})

test_that("metrics the null holds fixed are excluded from automatic selection", {
  # Nodes, Edges and Density have zero variance under a degree-preserving model,
  # so a z-score for them is undefined. Dropping them beats reporting NA rows
  # that look like failures.
  res <- metric_significance(g, n = 10, seed = 1, plot = FALSE)
  expect_false(any(c("Nodes", "Edges", "Density") %in% res$result$metric))
  expect_gt(nrow(res$result), 0)
})

test_that("a fixed metric requested explicitly is reported as NA, not dropped", {
  res <- metric_significance(g, metrics = c("Nodes", "Modularity"), n = 10,
                             seed = 1, plot = FALSE)
  expect_true("Nodes" %in% res$result$metric)
  expect_true(is.na(res$result$z[res$result$metric == "Nodes"]))
  expect_false(is.na(res$result$z[res$result$metric == "Modularity"]))
})

test_that("an unknown metric is rejected with the available names", {
  expect_error(metric_significance(g, metrics = "Smallworldness", n = 5,
                                   plot = FALSE),
               "Unknown or non-numeric metric")
})

test_that("a caller-supplied ensemble is used instead of a fresh one", {
  nulls <- null_model(g, model = "configuration", n = 12, seed = 5)
  res <- metric_significance(g, null = nulls, metrics = "Modularity",
                             plot = FALSE)
  expect_match(res$method, "12 null graphs")
  expect_match(res$method, "configuration")

  expect_error(metric_significance(g, null = list(), metrics = "Modularity"),
               "non-empty list")
})

test_that("the diagnostic plot builds and marks the observed value", {
  res <- metric_significance(g, metrics = c("Modularity",
                                            "Clustering_coefficient"),
                             n = 15, seed = 1, plot = TRUE)
  expect_s3_class(res$plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(res$plot))
})

# --- small_worldness() ------------------------------------------------------

test_that("small_worldness() separates a small-world graph from a random one", {
  # sigma >> 1 is the definition of small-world; a random graph is its own null,
  # so sigma should sit near 1. Asserted as an ordering plus a loose bound rather
  # than as a recorded number, since sigma depends on the ensemble.
  set.seed(21)
  sw <- igraph::sample_smallworld(1, 80, 4, 0.05)
  er <- igraph::sample_gnp(80, igraph::edge_density(sw), directed = FALSE)

  sw_res <- small_worldness(sw, n = 10, seed = 1)$result
  er_res <- small_worldness(er, n = 10, seed = 1)$result

  expect_gt(sw_res$sigma, 2)
  expect_lt(er_res$sigma, 2)
  expect_gt(sw_res$sigma, er_res$sigma)
})

test_that("small_worldness() reports its components and shape", {
  set.seed(22)
  sw <- igraph::sample_smallworld(1, 60, 4, 0.05)
  res <- small_worldness(sw, n = 8, seed = 1)

  expect_true(all(c("result", "graph", "method") %in% names(res)))
  expect_named(res$result, c("sigma", "C", "C_rand", "L", "L_rand", "n_null"))
  expect_equal(res$result$n_null, 8)
  # sigma must be reconstructible from the parts that were reported.
  expect_equal(res$result$sigma,
               (res$result$C / res$result$C_rand) /
                 (res$result$L / res$result$L_rand),
               tolerance = 1e-10)
})

test_that("small_worldness() returns NA rather than dividing by zero", {
  # A triangle-free graph has C = 0 and so does its degree-preserving null,
  # making C/C_rand a 0/0. NA is the honest answer.
  res <- small_worldness(ring_graph(), n = 5, seed = 1)
  expect_true(is.na(res$result$sigma))
})

# --- The degree-matched diffusion null --------------------------------------

test_that("permuted seeds match the real seeds' degrees where candidates exist", {
  set.seed(31)
  pa <- igraph::sample_pa(200, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(200))
  deg <- igraph::degree(pa)
  names(deg) <- igraph::V(pa)$name

  set.seed(9)
  seeds <- sample(names(deg), 10)
  sampler <- suppressWarnings(
    make_degree_matched_sampler(pa, seeds, setdiff(names(deg), seeds), 10)
  )

  draws <- replicate(200, deg[sampler()])
  per_seed_mean <- rowMeans(draws)

  # The match is per seed, so the right statistic is the per-seed discrepancy.
  # Most seeds of a scale-free graph have common degrees and are matched exactly;
  # the median absolute difference is therefore zero. The mean is dragged up by
  # the one hub in the seed set that has no comparable partner, which is the
  # documented limitation rather than a sampling fault.
  expect_equal(stats::median(abs(per_seed_mean - deg[seeds])), 0,
               tolerance = 0.5)

  # And the matched null must be much closer to the real seeds than a uniform
  # draw, which is the entire point of the correction.
  uniform_mean <- mean(replicate(
    200, mean(deg[sample(setdiff(names(deg), seeds), length(seeds))])
  ))
  expect_lt(abs(mean(per_seed_mean) - mean(deg[seeds])),
            abs(uniform_mean - mean(deg[seeds])))
})

test_that("match_pool = Inf reduces exactly to the uniform null", {
  # A structural equivalence, not a coincidence: with every candidate eligible,
  # each seed is drawn uniformly.
  g2 <- test_graph_gnp(n = 40, p = 0.1, seed = 41)
  seeds <- c("m1", "m2", "m3")

  set.seed(5)
  a <- network_diffusion_with_pvalues(g2, seeds, method = "rwr",
                                      n_permutations = 30, seed = 2,
                                      verbose = FALSE, null = "degree_matched",
                                      match_pool = Inf)
  set.seed(5)
  b <- network_diffusion_with_pvalues(g2, seeds, method = "rwr",
                                      n_permutations = 30, seed = 2,
                                      verbose = FALSE, null = "uniform")
  expect_equal(a, b)
})

test_that("a tighter match_pool tracks seed degree more closely", {
  set.seed(32)
  pa <- igraph::sample_pa(200, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(200))
  deg <- igraph::degree(pa)
  names(deg) <- igraph::V(pa)$name
  set.seed(9)
  seeds <- sample(names(deg), 10)

  discrepancy <- vapply(c(3, 30), function(mp) {
    sampler <- suppressWarnings(
      make_degree_matched_sampler(pa, seeds, setdiff(names(deg), seeds), mp)
    )
    draws <- replicate(200, deg[sampler()])
    mean(abs(rowMeans(draws) - deg[seeds]))
  }, numeric(1))

  expect_lt(discrepancy[1], discrepancy[2])
})

test_that("the degree-matched null is the default and both nulls are valid", {
  g2 <- test_graph_gnp(n = 40, p = 0.1, seed = 42)
  seeds <- c("m1", "m2", "m3")

  default_res <- network_diffusion_with_pvalues(g2, seeds, method = "rwr",
                                                n_permutations = 30, seed = 1,
                                                verbose = FALSE)
  explicit <- network_diffusion_with_pvalues(g2, seeds, method = "rwr",
                                             n_permutations = 30, seed = 1,
                                             verbose = FALSE,
                                             null = "degree_matched")
  expect_equal(default_res, explicit)

  for (res in list(default_res, explicit)) {
    expect_named(res, c("node", "score", "p_empirical", "null_mean", "null_sd",
                        "z", "n_valid"))
    expect_equal(nrow(res), igraph::vcount(g2))
    # Per-vertex floor: n_valid excludes the permutations that seeded the vertex.
    expect_true(all(res$p_empirical >= 1 / (res$n_valid + 1) - 1e-12))
    expect_true(all(res$p_empirical <= 1))
  }
})

test_that("the two nulls disagree on a graph where degree is confounded", {
  # If they agreed, the correction would be doing nothing. High-degree seeds on a
  # scale-free graph are exactly the case the degree-matched null exists for.
  set.seed(33)
  pa <- igraph::sample_pa(150, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(150))
  hubs <- names(sort(igraph::degree(pa), decreasing = TRUE))[1:5]

  u <- network_diffusion_with_pvalues(pa, hubs, method = "rwr",
                                      n_permutations = 60, seed = 1,
                                      verbose = FALSE, null = "uniform")
  d <- suppressWarnings(
    network_diffusion_with_pvalues(pa, hubs, method = "rwr",
                                   n_permutations = 60, seed = 1,
                                   verbose = FALSE, null = "degree_matched")
  )

  u <- u[order(u$node), ]
  d <- d[order(d$node), ]
  expect_false(isTRUE(all.equal(u$p_empirical, d$p_empirical)))

  # This used to assert that the uniform null finds at least as many significant
  # nodes as the matched one. That is not a theorem, and it stopped holding once
  # each vertex began being tested only against the permutations that did not
  # seed it: the degree-matched sampler draws from a narrower support and so
  # draws individual vertices more often (measured at up to 10 of 60
  # permutations, against 5 for uniform), so the exclusion removes a larger
  # share of its nulls' upper tails and lowers those p-values further. Both
  # p-values remain valid; the count comparison simply has two effects pulling
  # on it and no longer measures what it was added to measure.
  #
  # What it was added to measure was the degenerate sampler -- a null assembled
  # from ten vertices, which has almost no variance and so makes nearly any
  # observed score extreme. That is now asserted directly, on the variance
  # itself, which is both the mechanism and a quantity the function reports.
  expect_false(any(d$null_sd == 0, na.rm = TRUE))
  expect_gt(stats::median(d$null_sd, na.rm = TRUE),
            0.5 * stats::median(u$null_sd, na.rm = TRUE))
  # The sampler's own support is guarded in "the sampler explores the candidate
  # pool rather than one corner of it", which is the other half of this.
})

test_that("an unmatchable seed set warns rather than claiming an exact match", {
  set.seed(34)
  pa <- igraph::sample_pa(200, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(200))
  hubs <- names(sort(igraph::degree(pa), decreasing = TRUE))[1:6]

  # match_pool = NULL, so the pool-size warning does not fire and this asserts
  # the one it is about: the graph has no vertices comparable to these seeds.
  expect_warning(
    make_degree_matched_sampler(pa, hubs, setdiff(igraph::V(pa)$name, hubs),
                                NULL),
    "only approximately"
  )
})

test_that("a seed set too large for a permutation null errors clearly", {
  small <- test_graph_gnp(n = 10, p = 0.3, seed = 43)
  expect_error(
    network_diffusion_with_pvalues(small, igraph::V(small)$name[1:9],
                                   n_permutations = 5, verbose = FALSE),
    "too large a share of the graph"
  )
})

test_that("an invalid match_pool is rejected", {
  g2 <- test_graph_gnp(n = 40, p = 0.1, seed = 44)
  expect_error(
    network_diffusion_with_pvalues(g2, c("m1", "m2"), n_permutations = 5,
                                   verbose = FALSE, match_pool = 0),
    "at least 1"
  )
})

# --- Input contract ---------------------------------------------------------

test_that("the new functions honour the shared graph-input contract", {
  for (fn in list(null_model, metric_significance, small_worldness)) {
    expect_error(fn("not a graph"),
                 "Input 'graph' must be either an igraph object or a data.frame")
  }
  expect_length(null_model(as_edge_df(g), n = 2, seed = 1), 2)
})

# --- Weights under a name other than `weight` --------------------------------
#
# null_model() read only the literal attribute `weight`, so a graph annotated
# with any other name -- which is the normal case -- produced unweighted nulls
# while the caller had explicitly asked for weighted ones.

test_that("null_model() permutes a weight attribute under any name", {
  gw <- test_graph()
  set.seed(51)
  igraph::E(gw)$corr <- runif(igraph::ecount(gw), 0.1, 2)

  nulls <- null_model(gw, model = "rewire", n = 5, seed = 1, weights = "corr")
  for (null_g in nulls) {
    # Written under `weight` whatever the source was called, so consumers have
    # one name to ask for.
    expect_equal(igraph::edge_attr_names(null_g), "weight")
    expect_equal(sort(igraph::E(null_g)$weight), sort(igraph::E(gw)$corr))
  }
})

test_that("null_model() accepts a numeric weight vector", {
  gw <- test_graph()
  set.seed(52)
  w <- runif(igraph::ecount(gw), 0.1, 2)

  nulls <- null_model(gw, n = 3, seed = 1, weights = w)
  expect_equal(sort(igraph::E(nulls[[1]])$weight), sort(w))
})

test_that("rewired nulls carry no edge attribute left over from the observed graph", {
  # rewire() starts from the observed graph, so every observed edge attribute is
  # inherited -- now bound to edges that have since been moved. Those values are
  # not a permutation of anything and a weighted metric would read them as data.
  gw <- test_graph()
  set.seed(53)
  igraph::E(gw)$corr <- runif(igraph::ecount(gw), 0.1, 2)
  igraph::E(gw)$label <- "observed"

  unweighted <- null_model(gw, model = "rewire", n = 2, seed = 1,
                           shuffle_weights = FALSE)
  expect_length(igraph::edge_attr_names(unweighted[[1]]), 0)

  weighted <- null_model(gw, model = "rewire", n = 2, seed = 1,
                         weights = "corr")
  expect_equal(igraph::edge_attr_names(weighted[[1]]), "weight")
})

test_that("metric_significance() and small_worldness() accept non-standard weight names", {
  # Both used to tell null_model() only *whether* weights were wanted, never
  # which, so the nulls came back unweighted and the lookup for `weight` then
  # failed with an error naming an attribute the caller had never mentioned.
  gw <- test_graph()
  set.seed(54)
  igraph::E(gw)$corr <- runif(igraph::ecount(gw), 0.1, 2)

  res <- metric_significance(gw, weights = "corr", metrics = "Modularity",
                             n = 5, seed = 1, plot = FALSE)
  expect_equal(nrow(res$result), 1)
  expect_true(is.finite(res$result$observed))
  expect_match(res$method, "weighted by 'corr'")

  expect_no_error(
    metric_significance(gw, weights = igraph::E(gw)$corr,
                        metrics = "Modularity", n = 5, seed = 1, plot = FALSE)
  )
  expect_true(is.finite(small_worldness(gw, n = 5, seed = 1,
                                        weights = "corr")$result$C))
})

test_that("a caller-supplied unweighted ensemble is rejected, not silently compared", {
  gw <- test_graph()
  set.seed(55)
  igraph::E(gw)$corr <- runif(igraph::ecount(gw), 0.1, 2)

  bare <- null_model(gw, n = 3, seed = 1, shuffle_weights = FALSE)
  expect_error(
    metric_significance(gw, null = bare, weights = "corr",
                        metrics = "Modularity", plot = FALSE),
    "carry no 'weight' edge attribute"
  )
})

# --- The sampler's pool must be wide enough to be a null --------------------
#
# The per-seed nearest-degree pools this replaces matched degree tightly and
# were, as a null, nearly degenerate: ties broke by vertex index, so seeds of
# similar degree got near-identical pools and drew from them without
# replacement. 200 permutations of a 5-seed set on a 150-vertex graph touched
# 10 of 145 candidates. That was invisible until the permutation rows were
# aligned, because the p-values were scrambled anyway.

test_that("the sampler explores the candidate pool rather than one corner of it", {
  set.seed(71)
  pa <- igraph::sample_pa(150, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(150))
  hubs <- names(sort(igraph::degree(pa), decreasing = TRUE))[1:5]
  candidates <- setdiff(igraph::V(pa)$name, hubs)

  sampler <- suppressWarnings(
    make_degree_matched_sampler(pa, hubs, candidates, NULL)
  )
  draws <- replicate(200, sampler())

  # The old implementation reached 10 of 145. The exact figure depends on the
  # degree distribution, so this asserts the order of magnitude, not a number.
  expect_gt(length(unique(as.vector(draws))), 40)

  # A seed set is a set: no permutation may use one vertex twice, which would
  # quietly make the null's effective seed count smaller than the real one.
  expect_equal(sum(apply(draws, 2, anyDuplicated)), 0)

  # And a real seed must never appear in its own null.
  expect_length(intersect(as.vector(draws), hubs), 0)
})

test_that("ties are broken at random, not by vertex index", {
  # order() is stable, so the previous implementation always returned the
  # lowest-indexed of a set of equally good candidates. On a graph where every
  # candidate is equally good that made the sampler fully deterministic.
  ring <- igraph::make_ring(60)
  igraph::V(ring)$name <- paste0("r", seq_len(60))
  seeds <- c("r1", "r2")
  candidates <- setdiff(igraph::V(ring)$name, seeds)

  sampler <- make_degree_matched_sampler(ring, seeds, candidates, NULL)
  draws <- replicate(100, sampler())

  # Every vertex of a ring has degree 2, so all candidates are interchangeable
  # and the draws must spread over them rather than collapsing onto the first.
  expect_gt(length(unique(as.vector(draws))), 30)
})

test_that("the drawn set stays in seed order, so seed_weights stay attached", {
  # The sampler processes seeds in a random order to avoid giving the first one
  # the widest choice, but network_diffusion_with_pvalues() reuses
  # `seed_weights` against the returned vector positionally. If the shuffle
  # leaked into the return value, each weight would land on the wrong
  # replacement -- and nothing about the output would show it.
  set.seed(72)
  pa <- igraph::sample_pa(200, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(200))
  deg <- igraph::degree(pa)
  names(deg) <- igraph::V(pa)$name

  # Seeds of sharply different degree, so a permuted position can be traced
  # back to the seed it replaced.
  seeds <- c(names(sort(deg, decreasing = TRUE))[1],
             names(sort(deg))[1],
             names(sort(deg, decreasing = TRUE))[2])
  sampler <- suppressWarnings(
    make_degree_matched_sampler(pa, seeds, setdiff(names(deg), seeds), 5)
  )

  draws <- replicate(200, sampler())
  expect_equal(nrow(draws), length(seeds))

  # Position 2 replaces the lowest-degree seed, so its mean degree must stay
  # below those of positions 1 and 3. Under a leaked shuffle all three
  # positions would converge on the same mean.
  position_mean <- rowMeans(matrix(deg[as.vector(draws)], nrow = length(seeds)))
  expect_lt(position_mean[2], position_mean[1])
  expect_lt(position_mean[2], position_mean[3])
})

test_that("the degree-matched null is the more conservative of the two", {
  # The property the whole correction exists to provide, and the one the broken
  # sampler inverted: because permuted seeds have degrees like the real ones,
  # the null scores are comparable to the real ones, so *fewer* nodes clear the
  # threshold than under a uniform null that seeds mostly leaves.
  set.seed(73)
  pa <- igraph::sample_pa(400, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(400))
  seeds <- names(sort(igraph::degree(pa), decreasing = TRUE))[50:60]

  u <- network_diffusion_with_pvalues(pa, seeds, method = "rwr",
                                      n_permutations = 100, seed = 1,
                                      verbose = FALSE, null = "uniform")
  d <- suppressWarnings(
    network_diffusion_with_pvalues(pa, seeds, method = "rwr",
                                   n_permutations = 100, seed = 1,
                                   verbose = FALSE, null = "degree_matched")
  )

  expect_lte(sum(d$p_empirical < 0.05), sum(u$p_empirical < 0.05))
})

test_that("a pool too small for the seed count is rejected loudly", {
  # The anti-conservative regime: each permutation would consume most of the
  # bin, so consecutive permutations repeat and the null stops varying.
  set.seed(74)
  pa <- igraph::sample_pa(400, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(400))
  deg <- igraph::degree(pa)
  names(deg) <- igraph::V(pa)$name

  # Mid-degree seeds, so the graph *can* match them and the only warning in
  # play is the one about pool size. The hubs of a scale-free graph would also
  # trip the "only approximately degree-matched" warning, which is a different
  # claim about a different problem.
  seeds <- names(sort(deg, decreasing = TRUE))[60:69]
  candidates <- setdiff(names(deg), seeds)

  expect_warning(make_degree_matched_sampler(pa, seeds, candidates, 5),
                 "anti-conservative")

  # The automatic choice scales with the seed count and so must not warn.
  expect_no_warning(make_degree_matched_sampler(pa, seeds, candidates, NULL))
})

# --- Summarizing the null ensemble ------------------------------------------
#
# metric_significance() and small_worldness() both recompute every metric for
# every null graph, which is the expensive half of each and is embarrassingly
# parallel. They now share summarize_null_ensemble(), which distributes the work
# with future.apply under whatever plan the caller has set.

test_that("the ensemble is summarized in input order, one row per graph", {
  # The rows are bound positionally and each metric column is then treated as a
  # null distribution, so a reordering would be entirely invisible in the
  # output. Graphs with distinct edge counts make the order observable.
  ens <- lapply(c(3L, 7L, 5L, 9L, 4L), function(k) {
    gg <- igraph::make_full_graph(5)
    igraph::subgraph_from_edges(gg, igraph::E(gg)[seq_len(k)],
                                delete.vertices = FALSE)
  })

  rows <- summarize_null_ensemble(ens, NULL, "strength")

  expect_equal(nrow(rows), length(ens))
  expect_equal(rows$Edges, vapply(ens, igraph::ecount, numeric(1)))
})

test_that("an ensemble summary is unchanged by the plan it ran under", {
  # The parallel path cannot be exercised under the default sequential plan, and
  # the one parallel defect this package has had -- a package-level
  # future::plan() call that leaked worker connections -- was invisible locally
  # and only caught by CI. Forked futures are unavailable on Windows, hence the
  # skip rather than a hard requirement.
  skip_on_cran()
  skip_if_not(future::supportsMulticore(),
              "forked futures are unavailable on this platform")

  g <- test_graph(n = 40)
  nulls <- null_model(g, n = 6, seed = 2)

  sequential <- metric_significance(g, null = nulls, plot = FALSE)$result
  sw_sequential <- small_worldness(g, n = 6, seed = 2)$result

  plan_before <- future::plan()
  conns_before <- nrow(showConnections(all = FALSE))
  on.exit(future::plan(plan_before), add = TRUE)
  # suppressWarnings: parallelly warns when asked for more workers than it
  # detects cores, which is a property of the machine, not of netkit.
  suppressWarnings(future::plan("multicore", workers = 2))

  parallel <- metric_significance(g, null = nulls, plot = FALSE)$result
  sw_parallel <- small_worldness(g, n = 6, seed = 2)$result

  # Identical, not merely similar: Modularity is seeded inside
  # summarize_graph_metrics() through with_local_seed(), so distributing the
  # ensemble must not hand each worker its own random stream.
  expect_equal(parallel, sequential)
  expect_equal(sw_parallel, sw_sequential)

  future::plan(plan_before)
  expect_equal(class(future::plan()), class(plan_before))
  expect_equal(nrow(showConnections(all = FALSE)), conns_before)
})

test_that("neither function changes the caller's future plan", {
  # The regression guard that applies on every platform, parallel or not.
  g <- test_graph(n = 30)
  plan_before <- class(future::plan())
  conns_before <- nrow(showConnections(all = FALSE))

  metric_significance(g, metrics = "Modularity", n = 4, seed = 1, plot = FALSE)
  small_worldness(g, n = 4, seed = 1)

  expect_equal(class(future::plan()), plan_before)
  expect_equal(nrow(showConnections(all = FALSE)), conns_before)
})
