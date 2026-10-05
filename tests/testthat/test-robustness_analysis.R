g <- test_graph()

test_that("robustness_analysis() works with its default arguments", {
  # Regression guard: `removal_strategy` was left as its length-3 default vector
  # and never passed through match.arg(), so the documented default call failed
  # with "Invalid `removal_strategy`".
  expect_no_error(
    without_progress(robustness_analysis(g, steps = 5, n_reps = 2, plot = FALSE, seed = 1))
  )
})

test_that("the default strategy resolves to 'random'", {
  set.seed(1)
  default  <- without_progress(
    robustness_analysis(g, steps = 5, n_reps = 2, plot = FALSE, seed = 1))
  set.seed(1)
  explicit <- without_progress(
    robustness_analysis(g, removal_strategy = "random", steps = 5,
                        n_reps = 2, plot = FALSE, seed = 1))
  expect_equal(default$summary, explicit$summary)
})

test_that("the returned structure is as documented", {
  res <- robustness_analysis(g, removal_strategy = "degree", steps = 5,
                             n_reps = 1, plot = FALSE, seed = 1)

  expect_named(res, c("plot", "result", "all_results", "summary", "auc"))
  expect_null(res$plot)   # plot = FALSE, but the element is still present
  expect_s3_class(res$summary, "data.frame")
  expect_named(res$summary, c("rep", "removed", "removed_frac",
                              "lcc_size", "efficiency", "n_components"))
  expect_named(res$auc, c("lcc_size", "efficiency", "n_components"))
  for (a in res$auc) expect_true(is.numeric(a) && length(a) == 1)
})

test_that("all documented removal strategies run", {
  for (s in c("random", "degree", "betweenness")) {
    res <- without_progress(
      robustness_analysis(g, removal_strategy = s, steps = 5,
                          n_reps = 2, plot = FALSE, seed = 1))
    expect_type(res, "list")
    expect_gt(nrow(res$summary), 0)
  }
})

test_that("a numeric vertex attribute can drive removal order", {
  # `score` is attached by the fixture; this is the documented "any numeric
  # vertex attribute as priority" path.
  res <- robustness_analysis(g, removal_strategy = "score", steps = 5,
                             n_reps = 1, plot = FALSE, seed = 1)
  expect_type(res, "list")
  expect_gt(nrow(res$summary), 0)
})

test_that("an unknown strategy is rejected", {
  expect_error(
    robustness_analysis(g, removal_strategy = "not_an_attribute",
                        steps = 5, plot = FALSE, seed = 1),
    "Invalid `removal_strategy`"
  )
})

test_that("removing nodes shrinks the largest connected component", {
  res <- robustness_analysis(g, removal_strategy = "degree", steps = 10,
                             n_reps = 1, plot = FALSE, seed = 1)
  s <- res$summary[order(res$summary$removed), ]

  expect_true(all(diff(s$removed) > 0))
  expect_true(all(s$removed_frac >= 0 & s$removed_frac <= 1))
  # Targeting high-degree nodes must not grow the giant component.
  expect_lte(utils::tail(s$lcc_size, 1), utils::head(s$lcc_size, 1))
})

test_that("targeted attack is at least as damaging as random failure", {
  # The classic Albert-Jeong-Barabasi result: scale-free networks are robust to
  # random failure but fragile to degree-targeted attack, so the area under the
  # LCC curve should be no larger for the targeted strategy.
  rand <- without_progress(
    robustness_analysis(g, removal_strategy = "random", steps = 20,
                        n_reps = 10, plot = FALSE, seed = 1))
  targ <- robustness_analysis(g, removal_strategy = "degree", steps = 20,
                              n_reps = 1, plot = FALSE, seed = 1)

  expect_lt(targ$auc$lcc_size, rand$auc$lcc_size)
})

test_that("a seed makes random removal reproducible", {
  args <- list(graph = g, removal_strategy = "random", steps = 5,
               n_reps = 3, plot = FALSE, seed = 123)
  expect_equal(without_progress(do.call(robustness_analysis, args))$summary,
               without_progress(do.call(robustness_analysis, args))$summary)
})

test_that("AUC is finite for every metric under default arguments", {
  # Regression guard. Global efficiency divides by n(n-1), which is 0 once a single
  # vertex remains. With the default steps = 50, step_size becomes 1 for any graph
  # of about 51 nodes or fewer, so the loop did reach one vertex and auc$efficiency
  # came back NaN from the documented default call.
  for (n in c(20, 50, 51)) {
    gg <- test_graph(n = n, seed = 2)
    res <- without_progress(
      robustness_analysis(gg, removal_strategy = "degree", plot = FALSE, seed = 1)
    )
    expect_false(any(is.nan(res$summary$efficiency)), info = paste("n =", n))
    for (m in names(res$auc)) {
      expect_true(is.finite(res$auc[[m]]), info = paste("n =", n, "metric", m))
    }
  }
})

test_that("a collapsed network reports zero efficiency rather than NaN", {
  gg <- test_graph(n = 30, seed = 4)
  res <- robustness_analysis(gg, removal_strategy = "degree", steps = 40,
                             plot = FALSE, seed = 1)

  # steps > vcount forces step_size = 1, so the final step leaves one vertex.
  expect_equal(min(res$summary$lcc_size), 1)
  expect_equal(utils::tail(res$summary$efficiency, 1), 0)
})

test_that("an all-zero metric does not produce a NaN AUC or NaN plot values", {
  # Normalizing by max() is a divide-by-zero when a metric is zero throughout, as
  # global efficiency is for a graph with no edges.
  ge <- igraph::make_empty_graph(20, directed = FALSE)
  igraph::V(ge)$name <- paste0("e", seq_len(20))

  res <- draw_quietly(
    robustness_analysis(ge, removal_strategy = "degree", steps = 5,
                        plot = TRUE, seed = 1)
  )

  expect_equal(res$auc$efficiency, 0)
  for (m in names(res$auc)) expect_true(is.finite(res$auc[[m]]), info = m)

  built <- ggplot2::ggplot_build(res$plot)
  ys <- unlist(lapply(built$data, function(d) d$y))
  expect_true(all(is.finite(ys)))
})

test_that("metrics can be requested selectively", {
  res <- robustness_analysis(g, removal_strategy = "degree",
                             metrics = "lcc_size", steps = 5,
                             n_reps = 1, plot = FALSE, seed = 1)
  expect_named(res$auc, "lcc_size")
})

# --- Regression: the RNG stream ---------------------------------------------
#
# robustness_analysis() used to call set.seed(seed) unconditionally. With the
# documented default seed = NULL that runs set.seed(NULL), which re-seeds the
# generator from the clock: a seeded script was not reproducible, and the
# caller's own RNG stream was destroyed as a side effect. Neither symptom is
# visible by reading the function or its output.

test_that("the default call does not touch the caller's RNG stream", {
  g <- test_graph(n = 30)

  set.seed(42)
  without_progress(robustness_analysis(g, steps = 5, n_reps = 2, plot = FALSE,
                                       metrics = "lcc_size"))
  after_first <- runif(3)

  set.seed(42)
  without_progress(robustness_analysis(g, steps = 5, n_reps = 2, plot = FALSE,
                                       metrics = "lcc_size"))
  after_second <- runif(3)

  expect_equal(after_first, after_second)
})

test_that("an outer seed makes the random strategy itself reproducible", {
  g <- test_graph(n = 30)

  set.seed(7)
  a <- without_progress(robustness_analysis(g, steps = 5, n_reps = 3, plot = FALSE,
                                            metrics = "lcc_size"))
  set.seed(7)
  b <- without_progress(robustness_analysis(g, steps = 5, n_reps = 3, plot = FALSE,
                                            metrics = "lcc_size"))

  expect_equal(a$all_results, b$all_results)
})

test_that("an explicit seed is still honoured", {
  g <- test_graph(n = 30)
  a <- without_progress(robustness_analysis(g, steps = 5, n_reps = 3, plot = FALSE,
                                            metrics = "lcc_size", seed = 99))
  b <- without_progress(robustness_analysis(g, steps = 5, n_reps = 3, plot = FALSE,
                                            metrics = "lcc_size", seed = 99))
  expect_equal(a$all_results, b$all_results)
})

test_that("`result` aliases `summary`", {
  res <- robustness_analysis(test_graph(n = 30), removal_strategy = "degree",
                             steps = 5, plot = FALSE)
  expect_true(all(c("result", "summary") %in% names(res)))
  expect_identical(res$result, res$summary)
})
