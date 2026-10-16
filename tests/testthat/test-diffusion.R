# The diffusion cluster is netkit's most intricate code. The
# adjacency-normalize -> Laplacian -> Cholesky/transition-matrix block used to be
# implemented twice, once in network_diffusion() and once in prepare_diffusion();
# network_diffusion() now delegates to prepare_diffusion(), so there is a single
# code path. The numeric-pin test below is what guards the maths now that the two
# routes can no longer disagree with each other.

g     <- test_graph()
seeds <- c("n1", "n2", "n3")
methods <- c("laplacian", "heat", "rwr")

test_that("diffusion returns one score per node, sorted descending", {
  for (m in methods) {
    res <- network_diffusion(g, seeds, method = m)
    expect_s3_class(res, "tbl_df")
    expect_named(res, c("node", "score"))
    expect_setequal(res$node, igraph::V(g)$name)
    expect_equal(nrow(res), igraph::vcount(g))
    expect_type(res$score, "double")
    expect_false(any(is.na(res$score)))
    expect_equal(res$score, sort(res$score, decreasing = TRUE),
                 info = paste("scores not sorted for method", m))
  }
})

test_that("diffusion scores match pinned reference values", {
  # A 6-node ring seeded at one vertex. These numbers were captured from the
  # implementation at the point network_diffusion() and prepare_diffusion() were
  # merged into one code path, and verified bit-identical across both routes
  # before and after that merge. They exist to catch an unintended change to the
  # diffusion maths itself, which no structural test would notice.
  g <- igraph::make_ring(6)
  igraph::V(g)$name <- letters[1:6]

  expected <- list(
    laplacian = c(0.6456263174, 0.1393781993, 0.0313535080,
                  0.0129102680, 0.0313535080, 0.1393781993),
    heat      = c(0.4657761538, 0.2080108694, 0.0509457439,
                  0.0163106196, 0.0509457439, 0.2080108694),
    rwr       = c(0.4239985981, 0.1771410030, 0.0821182562,
                  0.0574828834, 0.0821182562, 0.1771410030)
  )

  for (m in names(expected)) {
    res <- network_diffusion(g, "a", method = m)
    res <- res[order(res$node), ]
    expect_equal(res$score, expected[[m]], tolerance = 1e-8,
                 info = paste("diffusion maths changed for method", m))
  }
})

test_that("diffusion respects the symmetry of a ring", {
  # On a ring seeded at one vertex, scores must be symmetric about the seed.
  g <- igraph::make_ring(6)
  igraph::V(g)$name <- letters[1:6]

  for (m in methods) {
    res <- network_diffusion(g, "a", method = m)
    s <- stats::setNames(res$score, res$node)
    expect_equal(s[["b"]], s[["f"]], tolerance = 1e-10)
    expect_equal(s[["c"]], s[["e"]], tolerance = 1e-10)
    # The seed holds the maximum and the antipode the minimum.
    expect_equal(names(which.max(s)), "a")
    expect_equal(names(which.min(s)), "d")
  }
})

test_that("normalizing a directed graph warns from either entry point", {
  # network_diffusion() used to carry this warning in its own inline kernel
  # block; now that it delegates, the warning lives in prepare_diffusion() and
  # must still surface through both.
  gd <- igraph::make_ring(6, directed = TRUE)
  igraph::V(gd)$name <- letters[1:6]

  expect_warning(network_diffusion(gd, "a", method = "rwr", normalize = TRUE),
                 "Normalization for directed graphs is not supported")
  expect_warning(prepare_diffusion(gd, method = "rwr", normalize = TRUE),
                 "Normalization for directed graphs is not supported")

  # No warning when normalization is not requested.
  expect_no_warning(prepare_diffusion(gd, method = "rwr", normalize = FALSE))
})

test_that("precompute() reproduces the un-precomputed result exactly", {
  # Both routes now share one kernel implementation, so this no longer detects
  # drift between two copies; it guards that `precompute =` is actually wired
  # through and honoured rather than silently ignored.
  for (m in methods) {
    direct <- network_diffusion(g, seeds, method = m)
    kernel <- prepare_diffusion(g, method = m)
    cached <- network_diffusion(g, seeds, method = m, precompute = kernel)

    direct <- direct[order(direct$node), ]
    cached <- cached[order(cached$node), ]

    expect_equal(cached$node, direct$node)
    expect_equal(cached$score, direct$score,
                 info = paste("precompute diverged for method", m))
  }
})

test_that("a precomputed kernel for the wrong method is refused", {
  kernel <- prepare_diffusion(g, method = "rwr")
  expect_error(
    network_diffusion(g, seeds, method = "laplacian", precompute = kernel),
    "do not match pre-computed data"
  )
})

test_that("prepare_diffusion() returns the documented kernel structure", {
  lap <- prepare_diffusion(g, method = "laplacian")
  expect_named(lap, c("L", "ch", "P", "use_sparse_P", "method", "weights_key"))
  # The key is what lets network_diffusion() reject a kernel built from
  # different weights instead of silently diffusing over the wrong matrix.
  expect_equal(lap$weights_key, "unweighted")
  expect_equal(lap$method, "laplacian")
  expect_false(is.null(lap$ch))   # Cholesky factor only for the Laplacian path
  expect_null(lap$P)
  expect_false(lap$use_sparse_P)

  rwr <- prepare_diffusion(g, method = "rwr")
  expect_null(rwr$ch)
  expect_false(is.null(rwr$P))    # transition matrix only for the RWR path
  expect_true(rwr$use_sparse_P)

  # The Laplacian is square, symmetric for an undirected graph, and sparse.
  expect_equal(dim(lap$L), c(igraph::vcount(g), igraph::vcount(g)))
  expect_s4_class(lap$L, "Matrix")
  expect_equal(Matrix::rowSums(lap$L), rep(0, igraph::vcount(g)),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("seed nodes score above non-seeds", {
  for (m in methods) {
    res <- network_diffusion(g, seeds, method = m)
    seed_mean    <- mean(res$score[res$node %in% seeds])
    nonseed_mean <- mean(res$score[!res$node %in% seeds])
    expect_gt(seed_mean, nonseed_mean)
  }
})

test_that("an invalid method is rejected", {
  expect_error(network_diffusion(g, seeds, method = "bogus"), "should be one of")
  expect_error(prepare_diffusion(g, method = "bogus"), "should be one of")
})

test_that("diffusion is deterministic across repeated calls", {
  expect_equal(network_diffusion(g, seeds, method = "laplacian"),
               network_diffusion(g, seeds, method = "laplacian"))
})

test_that("permutation p-values come back in the documented shape", {
  res <- network_diffusion_with_pvalues(
    g, seeds, method = "laplacian",
    n_permutations = 25, seed = 42, verbose = FALSE
  )

  # The exact column set: a stray `stringsAsFactors` column used to leak in
  # here, because tibble() has no such argument.
  expect_named(res, c("node", "score", "p_empirical", "null_mean", "null_sd",
                      "z", "n_valid"))
  expect_equal(nrow(res), igraph::vcount(g))
  expect_setequal(res$node, igraph::V(g)$name)
  expect_true(all(res$p_empirical >= 0 & res$p_empirical <= 1))
  expect_false(any(is.na(res$p_empirical)))
  expect_equal(res$p_empirical, sort(res$p_empirical))

  # Each vertex is tested only against the permutations that did not seed it,
  # so the count it was compared against is reported and is never larger than
  # the number of permutations run.
  expect_type(res$n_valid, "integer")
  expect_true(all(res$n_valid <= 25))
  expect_true(all(res$n_valid > 0))
  # The p-value floor is per-vertex for that reason.
  expect_true(all(res$p_empirical >= 1 / (res$n_valid + 1) - 1e-12))
})

test_that("the caller's future plan is left untouched and no connections leak", {
  # Regression guard. network_diffusion_with_pvalues() used to call
  # future::plan("multisession") internally, which (a) silently overrode the
  # user's global plan and (b) started worker processes that were never shut
  # down. The open worker sockets made R CMD check fail with
  # "checking examples ... ERROR / connections left open". It was invisible on a
  # sandboxed dev machine that cannot spawn subprocesses, and only showed up on CI.
  plan_before <- class(future::plan())
  conns_before <- nrow(showConnections(all = FALSE))

  network_diffusion_with_pvalues(
    g, seeds, method = "laplacian",
    n_permutations = 10, seed = 1, verbose = FALSE
  )

  expect_equal(class(future::plan()), plan_before)
  expect_equal(nrow(showConnections(all = FALSE)), conns_before)
})

test_that("p-values are reproducible for a fixed seed", {
  args <- list(graph = g, seed_nodes = seeds, method = "laplacian",
               n_permutations = 25, seed = 7, verbose = FALSE)
  expect_equal(do.call(network_diffusion_with_pvalues, args),
               do.call(network_diffusion_with_pvalues, args))
})

# --- Permutation null alignment ---------------------------------------------
#
# network_diffusion() returns its nodes sorted by descending score, so the real
# scores and every permutation come back in a different vertex order. Combining
# them by position gave each node the distribution of the k-th *ranked* null
# score instead of its own. The p-values stayed in range and the table still
# looked well-formed, which is why this survived: nothing about the output
# announced that the rows had been shuffled.

test_that("each node's null distribution belongs to that node", {
  # The direct test of the defect: reconstruct the null by hand, aligning by
  # name, and require the function to agree. Nothing here depends on the
  # diffusion maths -- only on which column ends up next to which node.
  set.seed(301)
  g <- igraph::sample_pa(80, power = 1.2, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(80))
  all_nodes <- igraph::V(g)$name
  seeds <- all_nodes[1:4]

  res <- network_diffusion_with_pvalues(g, seeds, method = "laplacian",
                                        n_permutations = 40, seed = 7,
                                        verbose = FALSE, null = "uniform")

  # Rebuild it independently: same kernel, same seeds, explicit name alignment.
  precomp <- prepare_diffusion(g, method = "laplacian", alpha = 0.7,
                               normalize = TRUE)
  real <- network_diffusion(g, seeds, precompute = precomp)
  real <- stats::setNames(real$score, real$node)[all_nodes]

  # The permuted seed sets must be drawn the way the function draws them.
  # future_lapply(future.seed = ) uses parallel-safe L'Ecuyer streams, not the
  # ambient RNG, so set.seed() + sample() would produce different sets and this
  # comparison would fail for a reason that has nothing to do with alignment.
  cand <- setdiff(all_nodes, seeds)
  perm_seed_sets <- future.apply::future_lapply(
    seq_len(40), function(i) sample(cand, 4), future.seed = 7
  )
  nulls <- vapply(perm_seed_sets, function(ps) {
    d <- network_diffusion(g, ps, precompute = precomp)
    stats::setNames(d$score, d$node)[all_nodes]
  }, numeric(length(all_nodes)))

  # A vertex is compared only against the permutations that did not seed it, so
  # the mask has to be reconstructed here too -- from the same drawn sets, but
  # independently of how the function builds it.
  seeded <- vapply(perm_seed_sets, function(ps) all_nodes %in% ps,
                   logical(length(all_nodes)))
  expected <- vapply(seq_along(all_nodes), function(i) {
    keep <- !seeded[i, ]
    (sum(nulls[i, keep] >= real[[i]]) + 1) / (sum(keep) + 1)
  }, numeric(1))

  got <- res$p_empirical[match(all_nodes, res$node)]
  expect_equal(unname(got), unname(expected))

  # The same reconstruction for the reported null summaries.
  expect_equal(res$n_valid[match(all_nodes, res$node)],
               as.integer(rowSums(!seeded)))
  exp_mean <- vapply(seq_along(all_nodes),
                     function(i) mean(nulls[i, !seeded[i, ]]), numeric(1))
  expect_equal(res$null_mean[match(all_nodes, res$node)], unname(exp_mean))
  exp_z <- (unname(real) - exp_mean) /
    vapply(seq_along(all_nodes),
           function(i) stats::sd(nulls[i, !seeded[i, ]]), numeric(1))
  expect_equal(res$z[match(all_nodes, res$node)], unname(exp_z))
})

test_that("p-values track the scores they test rather than piling at the ends", {
  # The symptom the misalignment produced, asserted as a property so it cannot
  # come back under a different mechanism. A node scoring well above its own
  # null must get a small p-value, so rank correlation is strongly negative;
  # the broken version scrambled the pairing and sat near -0.5, with most nodes
  # pinned at exactly the floor or exactly 1.
  set.seed(302)
  g <- igraph::sample_pa(150, power = 1.2, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(150))
  seeds <- igraph::V(g)$name[order(-igraph::degree(g))][1:6]

  res <- network_diffusion_with_pvalues(g, seeds, method = "laplacian",
                                        n_permutations = 100, seed = 11,
                                        verbose = FALSE, null = "uniform")

  expect_lt(stats::cor(res$score, res$p_empirical, method = "spearman"), -0.8)

  # Not a bimodal pile-up: fewer than half the nodes may sit at either extreme.
  at_floor <- sum(res$p_empirical <= 1 / 101 + 1e-12)
  at_one <- sum(res$p_empirical == 1)
  expect_lt(at_floor + at_one, 0.5 * nrow(res))
})

test_that("a seed node outranks a random node under its own null", {
  # A sanity check with an answer known in advance: the seeds are where the
  # signal starts, so they must be among the most significant nodes. Under the
  # misaligned version they frequently were not.
  set.seed(303)
  g <- igraph::sample_pa(100, power = 1.2, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(100))
  seeds <- igraph::V(g)$name[c(10, 20, 30)]

  res <- network_diffusion_with_pvalues(g, seeds, method = "rwr",
                                        n_permutations = 80, seed = 13,
                                        verbose = FALSE, null = "uniform")

  seed_p <- res$p_empirical[match(seeds, res$node)]
  expect_true(all(seed_p <= stats::median(res$p_empirical)))
})

# --- Batched multi-RHS diffusion ---------------------------------------------
#
# The permutation null solves its diffusions in batches -- one sparse solve for
# many seed vectors -- rather than one diffusion per permutation. The batching is
# an optimization only: it must not change a single number. These tests pin that,
# because a batch whose columns or rows had been reordered would produce output
# just as well-formed as a correct one.

test_that("diffusing a batch equals diffusing each column alone", {
  # The property that makes batching legitimate. Asserted for all three methods,
  # since each batches by a different mechanism: a multi-RHS Cholesky solve, a
  # Taylor expansion over a matrix, and a matrix-valued fixed-point iteration.
  set.seed(401)
  g <- igraph::sample_pa(70, power = 1.2, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(70))
  all_nodes <- igraph::V(g)$name

  seed_sets <- lapply(1:9, function(i) sample(all_nodes, 4))
  F0 <- matrix(
    vapply(seed_sets, function(ps) make_seed_vector(all_nodes, ps, NULL),
           numeric(70)),
    nrow = 70
  )

  for (m in c("laplacian", "heat")) {
    kern <- prepare_diffusion(g, method = m)

    batched <- as.matrix(diffuse_columns(F0, kern, method = m))
    one_at_a_time <- vapply(seq_len(ncol(F0)), function(j) {
      as.numeric(diffuse_columns(F0[, j, drop = FALSE], kern, method = m))
    }, numeric(70))

    # Bit-identical, not merely close: a multi-RHS Cholesky solve and a Taylor
    # expansion over a matrix are the same arithmetic as the column-by-column
    # versions, so any difference at all would mean something else changed.
    expect_equal(batched, one_at_a_time, ignore_attr = TRUE, tolerance = 0,
                 info = paste("batched diffusion differs for method", m))
  }
})

test_that("a batched random walk is at least as converged as a single-column one", {
  # RWR is the one method where batching is not identical arithmetic: the
  # iteration stops when the *slowest* column has converged, so every other
  # column is carried past its own stopping point. The claim that makes this
  # harmless is that extra iterations of a contraction can only move a column
  # closer to its fixed point -- asserted here directly, per column, rather than
  # assumed. The practical consequence is that `method = "rwr"` p-values are
  # reproducible to the convergence tolerance rather than to the bit, while
  # "laplacian" and "heat" are reproducible to the bit.
  set.seed(401)
  g <- igraph::sample_pa(70, power = 1.2, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(70))
  all_nodes <- igraph::V(g)$name

  F0 <- matrix(
    vapply(lapply(1:9, function(i) sample(all_nodes, 4)),
           function(ps) make_seed_vector(all_nodes, ps, NULL), numeric(70)),
    nrow = 70
  )

  kern <- prepare_diffusion(g, method = "rwr")
  restart <- 0.3
  residual <- function(f, f0) {
    sum(abs((1 - restart) * as.numeric(kern$P %*% f) + restart * f0 - f))
  }

  batched <- as.matrix(diffuse_columns(F0, kern, method = "rwr"))
  tol <- 1e-6 * max(1, max(Matrix::colSums(abs(F0))))

  for (j in seq_len(ncol(F0))) {
    single <- as.numeric(diffuse_columns(F0[, j, drop = FALSE], kern,
                                         method = "rwr"))
    expect_lte(residual(batched[, j], F0[, j]), residual(single, F0[, j]))
    # And the two agree to within the convergence tolerance, so the choice is
    # invisible at the precision the method claims.
    expect_lt(max(abs(batched[, j] - single)), tol)
  }
})

test_that("the batched result is invariant to how it is split into batches", {
  # Pins the batch boundaries specifically. A batch size that divides the column
  # count evenly and one that leaves a short final batch must agree, or the
  # recombination is dropping or duplicating columns.
  set.seed(402)
  g <- igraph::sample_gnp(50, 0.1, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(50))
  all_nodes <- igraph::V(g)$name
  kern <- prepare_diffusion(g, method = "laplacian")

  F0 <- matrix(
    vapply(lapply(1:10, function(i) sample(all_nodes, 3)),
           function(ps) make_seed_vector(all_nodes, ps, NULL), numeric(50)),
    nrow = 50
  )

  whole <- as.matrix(diffuse_columns(F0, kern, method = "laplacian"))
  for (size in c(1, 3, 5, 10)) {
    idx <- split(seq_len(10), ceiling(seq_len(10) / size))
    chunked <- do.call(cbind, lapply(idx, function(j) {
      as.matrix(diffuse_columns(F0[, j, drop = FALSE], kern,
                                method = "laplacian"))
    }))
    expect_equal(whole, chunked, ignore_attr = TRUE,
                 info = paste("batch size", size))
  }
})

test_that("diffuse_columns() is indifferent to the matrix representation", {
  # The right-hand side is coerced to a dgeMatrix because a base matrix takes a
  # much slower path through Matrix's sparse product. That is a speed choice, so
  # it has to be numerically free.
  set.seed(403)
  g <- igraph::sample_gnp(40, 0.12, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(40))
  F0 <- matrix(make_seed_vector(igraph::V(g)$name, c("v1", "v2"), NULL), ncol = 1)

  for (m in methods) {
    kern <- prepare_diffusion(g, method = m)
    from_base <- as.numeric(diffuse_columns(F0, kern, method = m))
    from_dge <- as.numeric(
      diffuse_columns(Matrix::Matrix(F0, sparse = FALSE), kern, method = m))
    expect_equal(from_base, from_dge, info = m)
  }
})

test_that("diffusion_batch_size() stays within its bounds", {
  # At least one column, never more than there are permutations, and narrowing
  # as the graph grows -- the memory bound the batch exists to respect.
  expect_equal(diffusion_batch_size(10, 3), 3L)
  expect_equal(diffusion_batch_size(1e6, 200), 1L)
  expect_lte(diffusion_batch_size(6000, 200), 200L)
  expect_gte(diffusion_batch_size(6000, 200), 1L)
  expect_lte(diffusion_batch_size(50000, 200), diffusion_batch_size(500, 200))
})

test_that("the diffusion kernel is built once per call", {
  # Observable proxy for a repeated factorization. prepare_diffusion() warns when
  # asked to normalize a directed graph, so the warning count is the number of
  # times the kernel was built. It used to be two: the real-score call ran before
  # the kernel existed and silently built one of its own, which for "laplacian"
  # meant a second sparse Cholesky factorization of the whole graph whose result
  # was then discarded.
  gd <- igraph::make_ring(12, directed = TRUE)
  igraph::V(gd)$name <- paste0("d", seq_len(12))

  seen <- character(0)
  withCallingHandlers(
    network_diffusion_with_pvalues(gd, c("d1", "d2"), method = "laplacian",
                                   n_permutations = 5, seed = 1,
                                   verbose = FALSE, null = "uniform"),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(sum(grepl("Normalization for directed graphs", seen)), 1L)
})

# --- Self-seeded permutations are excluded per node -------------------------
#
# A permutation that happened to seed vertex v deposits a large score at v, for
# the same reason the real seeds top the real run, and so says nothing about
# whether v is unusually well placed relative to the real seed set. Including it
# inflated the upper tail of v's own null and made its p-value conservative --
# worst for the vertices drawn most often, so not a uniform shift.

test_that("n_valid counts the permutations that did not seed the vertex", {
  set.seed(601)
  g <- igraph::sample_pa(60, power = 1.2, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(60))
  all_nodes <- igraph::V(g)$name
  seeds <- all_nodes[1:3]

  res <- network_diffusion_with_pvalues(g, seeds, method = "laplacian",
                                        n_permutations = 50, seed = 9,
                                        verbose = FALSE, null = "uniform")

  # Reconstruct the draws the way the function draws them, then count.
  cand <- setdiff(all_nodes, seeds)
  sets <- future.apply::future_lapply(seq_len(50), function(i) sample(cand, 3),
                                      future.seed = 9)
  times_seeded <- table(factor(unlist(sets), levels = all_nodes))

  expect_equal(res$n_valid[match(all_nodes, res$node)],
               as.integer(50 - times_seeded), ignore_attr = TRUE)

  # The real seeds are never drawn as permuted seeds, so nothing is excluded for
  # them. This is also why their p-values are not evidence of anything.
  expect_true(all(res$n_valid[res$node %in% seeds] == 50))
  # And at least one non-seed vertex must actually have been drawn, or the test
  # would pass without exercising the exclusion at all.
  expect_true(any(res$n_valid < 50))
})

test_that("excluding self-seeded permutations lowers p rather than raising it", {
  # Direction matters: the permutations removed are the ones in which the vertex
  # scored highest, so removing them can only make the observed score look more
  # extreme, never less. Asserted against the all-permutations p-value computed
  # by hand from the same run.
  set.seed(602)
  g <- igraph::sample_pa(80, power = 1.2, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(80))
  all_nodes <- igraph::V(g)$name
  seeds <- all_nodes[1:4]

  res <- network_diffusion_with_pvalues(g, seeds, method = "laplacian",
                                        n_permutations = 60, seed = 3,
                                        verbose = FALSE, null = "uniform")

  precomp <- prepare_diffusion(g, method = "laplacian")
  real <- network_diffusion(g, seeds, precompute = precomp)
  real <- stats::setNames(real$score, real$node)[all_nodes]
  cand <- setdiff(all_nodes, seeds)
  sets <- future.apply::future_lapply(seq_len(60), function(i) sample(cand, 4),
                                      future.seed = 3)
  nulls <- vapply(sets, function(ps) {
    d <- network_diffusion(g, ps, precompute = precomp)
    stats::setNames(d$score, d$node)[all_nodes]
  }, numeric(length(all_nodes)))

  p_including <- (rowSums(nulls >= real) + 1) / 61
  p_excluding <- res$p_empirical[match(all_nodes, res$node)]

  expect_true(all(p_excluding <= p_including + 1e-12))
  # And the change is real, not vacuous: some vertex must have moved.
  expect_true(any(p_excluding < p_including - 1e-12))
})

test_that("z ranks vertices that tie at the p-value floor", {
  # The p-value floor is 1 / (n_valid + 1), so with any realistic permutation
  # count the strongest vertices all report the same p and the ordering among
  # them would otherwise be arbitrary. z is returned to break that tie, and the
  # table is ordered by it.
  set.seed(603)
  g <- igraph::sample_pa(150, power = 1.2, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(150))
  seeds <- igraph::V(g)$name[order(-igraph::degree(g))][1:6]

  res <- suppressWarnings(
    network_diffusion_with_pvalues(g, seeds, method = "laplacian",
                                   n_permutations = 100, seed = 5,
                                   verbose = FALSE, null = "uniform"))

  expect_equal(res$p_empirical, sort(res$p_empirical))
  floor_rows <- res[res$p_empirical == min(res$p_empirical), ]
  expect_gt(nrow(floor_rows), 1)              # there is a tie to break
  expect_true(all(diff(floor_rows$z) <= 1e-12))  # ordered by descending z

  # z is the definition it claims to be, and tracks the score more closely than
  # a p-value pinned at its floor can.
  expect_equal(res$z, (res$score - res$null_mean) / res$null_sd)
  expect_gt(stats::cor(res$score, res$z, method = "spearman",
                       use = "complete.obs"),
            abs(stats::cor(res$score, res$p_empirical, method = "spearman")))
})

test_that("a vertex seeded in every permutation reports NA rather than a number", {
  # The degenerate case. With exactly as many candidates as seeds, every
  # permutation has to use all of them, so each candidate is seeded every time
  # and has no null left to be tested against. A p-value computed from nothing
  # would come out as 1 / (0 + 1) = 1, which looks like a result.
  g <- igraph::make_ring(6)
  igraph::V(g)$name <- letters[1:6]

  res <- network_diffusion_with_pvalues(
    g, seed_nodes = c("a", "b", "c"), method = "laplacian",
    n_permutations = 10, seed = 2, verbose = FALSE, null = "uniform")

  exhausted <- res[res$node %in% c("d", "e", "f"), ]
  expect_true(all(exhausted$n_valid == 0L))
  expect_true(all(is.na(exhausted$p_empirical)))
  expect_true(all(is.na(exhausted$z)))
  expect_true(all(is.na(exhausted$null_mean)))

  # The real seeds are never permuted, so their nulls are intact.
  real_seeds <- res[res$node %in% c("a", "b", "c"), ]
  expect_true(all(real_seeds$n_valid == 10L))
  expect_false(any(is.na(real_seeds$p_empirical)))

  # NA p-values sort last rather than being dropped or silently leading.
  expect_equal(nrow(res), 6L)
  expect_true(all(is.na(utils::tail(res$p_empirical, 3))))
})

# --- the degree-matched bin floor is derived, and the derivation holds -------
#
# match_pool = NULL used to default to 10 * n_seeds, a multiplier calibrated on
# five synthetic scale-free graphs. It is now derived from the self-seeding
# exclusion: a vertex seeded in a fraction r of permutations keeps n(1-r) of
# them, so its p-value floor rises by about 1/(1-r), and holding that within a
# tolerance t gives r <= t/(1+t) and hence a bin floor of n_seeds * (1+t)/t.
# At t = 0.1 that is 11 * n_seeds.
#
# These tests measure r by actually drawing permutations, rather than asserting
# the constant, so the derivation is checked against the sampler's behaviour.

seeding_rate <- function(graph, seeds, match_pool, n_draws = 2000) {
  sampler <- suppressWarnings(netkit:::make_degree_matched_sampler(
    graph, seeds, setdiff(igraph::V(graph)$name, seeds), match_pool))
  drawn <- unlist(replicate(n_draws, sampler(), simplify = FALSE))
  max(table(drawn)) / n_draws
}

test_that("the derived bin floor keeps the self-seeding rate within tolerance", {
  set.seed(5)
  g <- igraph::sample_pa(150, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(150))
  seeds <- paste0("v", c(2, 5, 9, 14, 21))

  # t = 0.1 bounds r at 0.1/1.1 = 0.0909. A little slack for sampling noise in
  # 2000 draws; the assertion is the bound, not a recorded value.
  derived <- 11 * length(seeds)
  expect_lt(seeding_rate(g, seeds, derived), 0.0909 + 0.02)
})

test_that("the self-seeding rate falls as the bin floor rises", {
  # Monotonicity is the substance of the derivation: r ~ n_seeds / B. Asserting
  # the trend rather than the values keeps this robust to the binning details.
  set.seed(5)
  g <- igraph::sample_pa(150, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(150))
  seeds <- paste0("v", c(2, 5, 9, 14, 21))

  rates <- vapply(c(5, 20, 55, 110),
                  function(mp) seeding_rate(g, seeds, mp), numeric(1))
  expect_true(all(diff(rates) < 0))
  # A pool at the seed-set size is the degenerate case the floor exists to
  # prevent: the same few vertices in every permutation.
  expect_gt(rates[[1]], 3 * rates[[4]])
})

test_that("a match_pool below the derived floor warns about the p-value floor", {
  set.seed(5)
  g <- igraph::sample_pa(150, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(150))
  seeds <- paste0("v", c(2, 5, 9, 14, 21, 28, 35, 42))

  # The warning now states the consequence -- how much of each vertex's own null
  # the exclusion removes -- rather than only that the number is small.
  expect_warning(
    netkit:::make_degree_matched_sampler(
      g, seeds, setdiff(igraph::V(g)$name, seeds), match_pool = 10),
    "p-value floor"
  )
})

test_that("match_pool = NULL derives the floor from the seed count", {
  # 11 * n_seeds, capped at the candidate count. Asserted through behaviour: the
  # derived pool must give a rate within tolerance for seed sets of different
  # sizes, which a fixed pool would not.
  set.seed(8)
  g <- igraph::sample_pa(400, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(400))

  for (k in c(3, 10)) {
    seeds <- paste0("v", seq(2, by = 7, length.out = k))
    auto <- suppressWarnings(netkit:::make_degree_matched_sampler(
      g, seeds, setdiff(igraph::V(g)$name, seeds), match_pool = NULL))
    drawn <- unlist(replicate(1500, auto(), simplify = FALSE))
    expect_lt(max(table(drawn)) / 1500, 0.0909 + 0.03)
  }
})
