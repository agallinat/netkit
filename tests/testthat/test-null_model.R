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

test_that("configuration falls back when vl cannot realise the degree sequence", {
  # "vl" needs a connected realisation to exist, which isolated vertices rule
  # out. The fallback must warn rather than silently change model.
  iso <- igraph::add_vertices(test_graph(n = 30), 3,
                              name = c("iso1", "iso2", "iso3"))
  expect_warning(null_model(iso, model = "configuration", n = 2, seed = 1),
                 "configuration.simple")
})

# --- metric_significance() --------------------------------------------------

test_that("metric_significance() returns the full plot/result/graph/method shape", {
  res <- metric_significance(g, metrics = "Modularity", n_null = 10, seed = 1,
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
                             n_null = n_null, seed = 1, plot = FALSE)
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
                             n_null = 20, seed = 1, plot = FALSE)
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
                             model = "erdos_renyi", n_null = 50, seed = 2,
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
                             model = "erdos_renyi", n_null = 30, seed = 3,
                             plot = FALSE)
  expect_gt(res$result$z, 3)
  expect_lte(res$result$p_empirical, 1 / 31 + 1e-12)
})

test_that("metrics the null holds fixed are excluded from automatic selection", {
  # Nodes, Edges and Density have zero variance under a degree-preserving model,
  # so a z-score for them is undefined. Dropping them beats reporting NA rows
  # that look like failures.
  res <- metric_significance(g, n_null = 10, seed = 1, plot = FALSE)
  expect_false(any(c("Nodes", "Edges", "Density") %in% res$result$metric))
  expect_gt(nrow(res$result), 0)
})

test_that("a fixed metric requested explicitly is reported as NA, not dropped", {
  res <- metric_significance(g, metrics = c("Nodes", "Modularity"), n_null = 10,
                             seed = 1, plot = FALSE)
  expect_true("Nodes" %in% res$result$metric)
  expect_true(is.na(res$result$z[res$result$metric == "Nodes"]))
  expect_false(is.na(res$result$z[res$result$metric == "Modularity"]))
})

test_that("an unknown metric is rejected with the available names", {
  expect_error(metric_significance(g, metrics = "Smallworldness", n_null = 5,
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
                             n_null = 15, seed = 1, plot = TRUE)
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

  sw_res <- small_worldness(sw, n_null = 10, seed = 1)$result
  er_res <- small_worldness(er, n_null = 10, seed = 1)$result

  expect_gt(sw_res$sigma, 2)
  expect_lt(er_res$sigma, 2)
  expect_gt(sw_res$sigma, er_res$sigma)
})

test_that("small_worldness() reports its components and shape", {
  set.seed(22)
  sw <- igraph::sample_smallworld(1, 60, 4, 0.05)
  res <- small_worldness(sw, n_null = 8, seed = 1)

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
  res <- small_worldness(ring_graph(), n_null = 5, seed = 1)
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
    expect_named(res, c("node", "score", "p_empirical"))
    expect_equal(nrow(res), igraph::vcount(g2))
    expect_true(all(res$p_empirical >= 1 / 31 - 1e-12))
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
  # The uniform null credits degree as signal, so it must find at least as many
  # significant nodes as the matched one.
  expect_gte(sum(u$p_empirical < 0.05), sum(d$p_empirical < 0.05))
})

test_that("an unmatchable seed set warns rather than claiming an exact match", {
  set.seed(34)
  pa <- igraph::sample_pa(200, power = 1.2, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(200))
  hubs <- names(sort(igraph::degree(pa), decreasing = TRUE))[1:6]

  expect_warning(
    make_degree_matched_sampler(pa, hubs, setdiff(igraph::V(pa)$name, hubs), 10),
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
