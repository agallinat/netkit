g <- test_graph()
targets <- c("n10", "n11", "n12")

test_that("greedy_seed_selection() returns the documented structure", {
  res <- without_leaf_warning(
    greedy_seed_selection(g, target_nodes = targets, k = 3,
                          method = "laplacian", plot = FALSE))

  expect_named(res, c("selected_seeds", "final_target_score",
                      "scores_at_each_step", "objective_at_each_step",
                      "result", "plot", "method"))
  expect_type(res$selected_seeds, "character")
  expect_length(res$selected_seeds, 3L)
  expect_length(res$scores_at_each_step, 3L)
  expect_true(is.numeric(res$final_target_score))
})

test_that("selected seeds are real graph vertices and exclude the targets", {
  res <- without_leaf_warning(
    greedy_seed_selection(g, target_nodes = targets, k = 3,
                          method = "laplacian", plot = FALSE))

  expect_true(all(res$selected_seeds %in% igraph::V(g)$name))
  expect_false(any(duplicated(res$selected_seeds)))
  # candidate_nodes defaults to setdiff(all_nodes, target_nodes).
  expect_length(intersect(res$selected_seeds, targets), 0L)
})

test_that("the greedy objective improves monotonically", {
  res <- without_leaf_warning(
    greedy_seed_selection(g, target_nodes = targets, k = 4,
                          method = "laplacian", plot = FALSE))

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
  implicit <- without_leaf_warning(
    greedy_seed_selection(g, target_nodes = targets, k = 2, plot = FALSE))
  explicit <- without_leaf_warning(
    greedy_seed_selection(g, target_nodes = targets, k = 2,
                          method = "laplacian", plot = FALSE))
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
    res <- without_leaf_warning(
      greedy_seed_selection(g, target_nodes = targets, k = 2,
                            method = m, plot = FALSE))
    expect_length(res$selected_seeds, 2L)
  }
})

# --- The linear objective makes the greedy search redundant ------------------
#
# Diffusion is a linear operator, so the target score of a seed set is the sum
# of its members' individual contributions. The default objective is therefore
# modular, every step ranks the candidates identically, and the result is the
# top-k single-seed contributions. These tests pin that as a fact about the
# method rather than leaving it an unstated assumption, because the whole design
# of `objective = "coverage"` exists to escape it.

test_that("the target score of a seed set is the sum of its members' scores", {
  # The property the entire contribution-matrix implementation rests on. Checked
  # against full diffusions of the combined set, not against the matrix.
  set.seed(501)
  gg <- igraph::sample_pa(80, power = 1.2, directed = FALSE)
  igraph::V(gg)$name <- paste0("v", seq_len(80))
  tg <- c("v3", "v17", "v40")
  members <- c("v50", "v60", "v70")

  for (m in c("laplacian", "heat")) {
    kern <- prepare_diffusion(gg, method = m)
    score <- function(S) {
      d <- network_diffusion(gg, S, method = m, precompute = kern)
      sum(d$score[d$node %in% tg])
    }
    expect_equal(score(members), sum(vapply(members, score, numeric(1))),
                 info = m)
  }

  # "rwr" is additive too -- it is the same linear algebra -- but it is solved
  # iteratively, so it is additive only to its convergence tolerance.
  kern <- prepare_diffusion(gg, method = "rwr")
  score_rwr <- function(S) {
    d <- network_diffusion(gg, S, method = "rwr", precompute = kern)
    sum(d$score[d$node %in% tg])
  }
  expect_equal(score_rwr(members),
               sum(vapply(members, score_rwr, numeric(1))),
               tolerance = 1e-6)
})

test_that("the 'total' objective selects exactly the top-k single contributions", {
  set.seed(502)
  gg <- igraph::sample_pa(100, power = 1.2, directed = FALSE)
  igraph::V(gg)$name <- paste0("v", seq_len(100))
  tg <- c("v5", "v25")
  cands <- setdiff(igraph::V(gg)$name, tg)
  kern <- prepare_diffusion(gg, method = "laplacian")

  singles <- vapply(cands, function(cd) {
    d <- network_diffusion(gg, cd, method = "laplacian", precompute = kern)
    sum(d$score[d$node %in% tg])
  }, numeric(1))

  res <- suppressWarnings(
    greedy_seed_selection(gg, tg, k = 5, method = "laplacian", plot = FALSE))

  expect_equal(res$selected_seeds,
               names(sort(singles, decreasing = TRUE))[1:5])
  # And the running total grows by each pick's own contribution, which is what
  # modularity means.
  expect_equal(res$scores_at_each_step,
               cumsum(sort(singles, decreasing = TRUE)[1:5]),
               ignore_attr = TRUE)
})

test_that("the 'total' objective reports its own redundancy", {
  res <- suppressWarnings(
    greedy_seed_selection(g, target_nodes = targets, k = 3, plot = FALSE))
  expect_match(res$method, "identical to the top-3")
  # Under the linear objective the two step vectors are the same quantity.
  expect_equal(res$objective_at_each_step, res$scores_at_each_step)
})

# --- The saturating objective -----------------------------------------------

test_that("'coverage' gains diminish and are bounded by the target count", {
  set.seed(503)
  gg <- igraph::sample_pa(150, power = 1.2, directed = FALSE)
  igraph::V(gg)$name <- paste0("v", seq_len(150))
  tg <- c("v4", "v30", "v77", "v120")

  res <- suppressWarnings(
    greedy_seed_selection(gg, tg, k = 6, plot = FALSE, objective = "coverage"))

  # Monotone, submodular, and bounded: one unit of credit per target at most.
  expect_true(all(diff(res$objective_at_each_step) >= 0))
  expect_lt(max(res$objective_at_each_step), length(tg))
  gains <- diff(c(0, res$objective_at_each_step))
  # Diminishing returns overall: the last pick must buy less than the first.
  expect_lt(gains[length(gains)], gains[1])
})

test_that("'coverage' spreads across targets while 'total' concentrates", {
  # The behavioural difference the objective exists to produce, asserted on the
  # per-target distribution rather than on which vertices get picked. A star of
  # targets each with its own private leaf: the total objective has no reason to
  # serve more than the cheapest targets, coverage has every reason to serve all.
  edges <- c(
    "hub", "t1", "hub", "t2", "hub", "t3",
    "t1", "a1", "t1", "a2", "t1", "a3",
    "t2", "b1", "t3", "c1"
  )
  gg <- igraph::make_graph(edges, directed = FALSE)
  tg <- c("t1", "t2", "t3")

  spread <- function(objective) {
    res <- suppressWarnings(
      greedy_seed_selection(gg, tg, k = 3, plot = FALSE, objective = objective))
    kern <- prepare_diffusion(gg, method = "laplacian")
    d <- network_diffusion(gg, res$selected_seeds, method = "laplacian",
                           precompute = kern)
    per_target <- d$score[match(tg, d$node)]
    # Coefficient of variation: low means the signal is shared out evenly.
    stats::sd(per_target) / mean(per_target)
  }

  expect_lt(spread("coverage"), spread("total"))
})

test_that("a coverage selection never beats a total selection on the total", {
  # The trade-off, stated as a property. "total" maximizes the summed score by
  # construction, so coverage can only match or lose on that scale -- and if it
  # ever won, the total objective would not be doing its job.
  set.seed(504)
  gg <- igraph::sample_pa(120, power = 1.2, directed = FALSE)
  igraph::V(gg)$name <- paste0("v", seq_len(120))
  tg <- c("v6", "v44", "v90")

  tot <- suppressWarnings(
    greedy_seed_selection(gg, tg, k = 4, plot = FALSE, objective = "total"))
  cov <- suppressWarnings(
    greedy_seed_selection(gg, tg, k = 4, plot = FALSE, objective = "coverage"))

  expect_lte(cov$final_target_score, tot$final_target_score + 1e-12)
})

test_that("an unknown objective is rejected eagerly", {
  expect_error(
    greedy_seed_selection(g, target_nodes = targets, k = 2, plot = FALSE,
                          objective = "bogus"),
    "should be one of"
  )
})

# --- The leaf-selection warning ---------------------------------------------

test_that("a selection dominated by degree-1 vertices warns", {
  # Symmetric normalization weights an edge into a degree-d vertex by
  # 1/sqrt(d), so a leaf hanging off a target is the most efficient injection
  # point in the graph. Correct, and usually not what the caller wanted.
  set.seed(505)
  gg <- igraph::sample_pa(200, power = 1.2, directed = FALSE)
  igraph::V(gg)$name <- paste0("v", seq_len(200))
  tg <- igraph::V(gg)$name[c(3, 17, 40)]

  expect_warning(greedy_seed_selection(gg, tg, k = 5, plot = FALSE),
                 "have degree 1")

  res <- suppressWarnings(greedy_seed_selection(gg, tg, k = 5, plot = FALSE))
  expect_gt(mean(igraph::degree(gg, v = res$selected_seeds) == 1), 0.5)
})

test_that("restricting candidates to non-leaves silences the warning", {
  # The remedy the warning names has to actually work.
  set.seed(505)
  gg <- igraph::sample_pa(200, power = 1.2, directed = FALSE)
  igraph::V(gg)$name <- paste0("v", seq_len(200))
  tg <- igraph::V(gg)$name[c(3, 17, 40)]
  interior <- setdiff(igraph::V(gg)$name[igraph::degree(gg) > 1], tg)

  expect_no_warning(
    greedy_seed_selection(gg, tg, candidate_nodes = interior, k = 5,
                          plot = FALSE)
  )
})

test_that("k is capped at the number of candidates", {
  res <- greedy_seed_selection(g, target_nodes = targets,
                               candidate_nodes = c("n20", "n21"), k = 5,
                               plot = FALSE)
  expect_length(res$selected_seeds, 2L)
})

test_that("an empty candidate set is an error, not an empty selection", {
  expect_error(
    greedy_seed_selection(g, target_nodes = targets,
                          candidate_nodes = character(0), k = 2, plot = FALSE),
    "No candidate nodes"
  )
})

# --- The per-seed diagnostic table ------------------------------------------

test_that("result has one row per selected seed with the documented columns", {
  res <- without_leaf_warning(
    greedy_seed_selection(g, target_nodes = targets, k = 4, plot = FALSE))

  expect_s3_class(res$result, "data.frame")
  expect_equal(nrow(res$result), 4L)
  expect_named(res$result, c("step", "node", "degree", "gain", "gain_z",
                             "cumulative_score", "objective"))
  expect_equal(res$result$node, res$selected_seeds)
  expect_equal(res$result$step, 1:4)
  expect_equal(res$result$degree,
               as.integer(igraph::degree(g, v = res$selected_seeds)),
               ignore_attr = TRUE)
  # The table and the step vectors are two views of one computation.
  expect_equal(res$result$cumulative_score, res$scores_at_each_step)
  expect_equal(res$result$objective, res$objective_at_each_step)
})

test_that("gain is the increment the pick bought", {
  # Under the linear objective the gain is the pick's own single-seed
  # contribution, so the cumulative score is the running sum of the gains. That
  # identity is the clearest statement of the modularity the help page warns
  # about.
  res <- without_leaf_warning(
    greedy_seed_selection(g, target_nodes = targets, k = 4, plot = FALSE))

  expect_equal(res$result$cumulative_score, cumsum(res$result$gain))
  # Greedy order: each gain is no larger than the one before it.
  expect_true(all(diff(res$result$gain) <= 1e-12))
})

test_that("gain_z standardizes against comparable degree, not against everything", {
  # The point of the column: a leaf's raw gain is large because it is a leaf. If
  # the z-score were taken over all candidates it would simply restate the gain,
  # so it must not correlate perfectly with it.
  set.seed(701)
  gg <- igraph::sample_pa(200, power = 1.2, directed = FALSE)
  igraph::V(gg)$name <- paste0("v", seq_len(200))
  tg <- igraph::V(gg)$name[c(3, 17, 40)]

  res <- suppressWarnings(
    greedy_seed_selection(gg, tg, k = 6, plot = FALSE))

  expect_equal(nrow(res$result), 6L)
  expect_false(isTRUE(all.equal(rank(res$result$gain),
                                rank(res$result$gain_z))))

  # Concretely: picks with identical gains get different z-scores, because what
  # they are compared against differs.
  same_gain <- res$result[abs(res$result$gain - res$result$gain[1]) < 1e-12, ]
  skip_if(nrow(same_gain) < 2, "fixture has no tied gains")
  expect_gt(stats::sd(same_gain$gain_z), 0)
})

test_that("a high-degree pick is not credited for its degree", {
  # Direct check of the separation. A star's centre delivers far more signal
  # than any of its leaves, but it is the only vertex of its degree, and
  # compared against the nearest-degree candidates -- the leaves -- it should
  # stand out. Compared on raw gain alone the two claims are indistinguishable.
  star <- igraph::make_star(30, mode = "undirected", center = 1)
  igraph::V(star)$name <- paste0("s", seq_len(30))

  res <- greedy_seed_selection(star, target_nodes = "s30", k = 2,
                               plot = FALSE)
  # The centre is adjacent to the target and is picked first.
  expect_equal(res$result$node[1], "s1")
  expect_gt(res$result$gain[1], res$result$gain[2])
})

test_that("gain_z is computed among the candidates still available", {
  # Documented behaviour, asserted so it is a choice rather than an accident:
  # the pool shrinks as picks are made, so identical gains earn rising z-scores
  # once the comparable alternatives have been consumed.
  #
  # Built by hand rather than sampled: six interchangeable leaves on one target
  # give exactly tied gains, which is the case the assertion needs and which no
  # random graph guarantees.
  g2 <- igraph::make_graph(
    c(c(rbind("t1", paste0("leaf", 1:6))),
      "t1", "h1", "h1", "h2", "h2", "h3", "h3", "t1"),
    directed = FALSE)

  res <- suppressWarnings(
    greedy_seed_selection(g2, target_nodes = "t1", k = 4, plot = FALSE))

  # The four picks are leaves, so their gains are identical to the bit.
  expect_equal(res$result$gain, rep(res$result$gain[1], 4))
  expect_true(all(res$result$degree == 1))
  # Same gain, strictly increasing z as the pool of equally good leaves empties.
  expect_true(all(diff(res$result$gain_z) > 0))
})

test_that("degree_matched_gain_z() pools by degree radius, not by rank", {
  # Unit-level: ties in degree distance must be taken whole. Resolving them by
  # index is what made the diffusion null's degree-matched sampler degenerate,
  # so the same mistake is excluded here by construction.
  gains <- c(10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 1)
  degrees <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 100, 100)

  # The first candidate is compared against the ten degree-5 candidates only;
  # the two degree-100 outliers are outside the radius that admits 10.
  z <- degree_matched_gain_z(gains, degrees, 1L, min_pool = 10L)
  pool <- gains[degrees == 5]
  expect_equal(z, (gains[1] - mean(pool)) / stats::sd(pool))

  # With no spread among comparable candidates there is nothing to standardize.
  expect_true(is.na(degree_matched_gain_z(rep(1, 12), degrees, 1L,
                                          min_pool = 10L)))
})
