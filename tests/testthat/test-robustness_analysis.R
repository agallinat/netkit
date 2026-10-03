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

  expect_named(res, c("plot", "all_results", "summary", "auc"))
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

test_that("metrics can be requested selectively", {
  res <- robustness_analysis(g, removal_strategy = "degree",
                             metrics = "lcc_size", steps = 5,
                             n_reps = 1, plot = FALSE, seed = 1)
  expect_named(res$auc, "lcc_size")
})
