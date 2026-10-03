g <- test_graph()
targets <- c("n10", "n11", "n12")

test_that("greedy_seed_selection() returns the documented structure", {
  res <- greedy_seed_selection(g, target_nodes = targets, k = 3,
                               method = "laplacian", plot = FALSE)

  expect_named(res, c("selected_seeds", "final_target_score",
                      "scores_at_each_step", "plot"))
  expect_type(res$selected_seeds, "character")
  expect_length(res$selected_seeds, 3L)
  expect_length(res$scores_at_each_step, 3L)
  expect_true(is.numeric(res$final_target_score))
})

test_that("selected seeds are real graph vertices and exclude the targets", {
  res <- greedy_seed_selection(g, target_nodes = targets, k = 3,
                               method = "laplacian", plot = FALSE)

  expect_true(all(res$selected_seeds %in% igraph::V(g)$name))
  expect_false(any(duplicated(res$selected_seeds)))
  # candidate_nodes defaults to setdiff(all_nodes, target_nodes).
  expect_length(intersect(res$selected_seeds, targets), 0L)
})

test_that("the greedy objective improves monotonically", {
  res <- greedy_seed_selection(g, target_nodes = targets, k = 4,
                               method = "laplacian", plot = FALSE)

  # Each greedy step adds the single best remaining seed, so the running target
  # score can never decrease.
  expect_true(all(diff(res$scores_at_each_step) >= 0))
  expect_equal(res$final_target_score,
               utils::tail(res$scores_at_each_step, 1),
               ignore_attr = TRUE)
})

test_that("the default method resolves to 'laplacian'", {
  # `method` was never passed through match.arg(), so the default silently
  # resolved downstream; it is now resolved eagerly.
  implicit <- greedy_seed_selection(g, target_nodes = targets, k = 2, plot = FALSE)
  explicit <- greedy_seed_selection(g, target_nodes = targets, k = 2,
                                    method = "laplacian", plot = FALSE)
  expect_equal(implicit$selected_seeds, explicit$selected_seeds)
})

test_that("an invalid method is rejected eagerly", {
  expect_error(
    greedy_seed_selection(g, target_nodes = targets, k = 2, method = "bogus",
                          plot = FALSE),
    "should be one of"
  )
})

test_that("candidate_nodes restricts the search space", {
  candidates <- c("n20", "n21", "n22")
  res <- greedy_seed_selection(g, target_nodes = targets,
                               candidate_nodes = candidates, k = 2,
                               method = "laplacian", plot = FALSE)
  expect_true(all(res$selected_seeds %in% candidates))
})

test_that("targets absent from the graph are an error", {
  expect_error(
    greedy_seed_selection(g, target_nodes = c("nope1", "nope2"), k = 2,
                          plot = FALSE),
    "No target nodes found"
  )
})

test_that("all diffusion methods are usable as the objective", {
  for (m in c("laplacian", "heat", "rwr")) {
    res <- greedy_seed_selection(g, target_nodes = targets, k = 2,
                                 method = m, plot = FALSE)
    expect_length(res$selected_seeds, 2L)
  }
})
