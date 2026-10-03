# The diffusion cluster is netkit's most intricate code. The
# adjacency-normalise -> Laplacian -> Cholesky/transition-matrix block used to be
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
  expect_named(lap, c("L", "ch", "P", "use_sparse_P", "method"))
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

  # Exactly three columns: a stray `stringsAsFactors` column used to leak in
  # here, because tibble() has no such argument.
  expect_named(res, c("node", "score", "p_empirical"))
  expect_equal(nrow(res), igraph::vcount(g))
  expect_setequal(res$node, igraph::V(g)$name)
  expect_true(all(res$p_empirical >= 0 & res$p_empirical <= 1))
  expect_false(any(is.na(res$p_empirical)))
  expect_equal(res$p_empirical, sort(res$p_empirical))
})

test_that("p-values are reproducible for a fixed seed", {
  args <- list(graph = g, seed_nodes = seeds, method = "laplacian",
               n_permutations = 25, seed = 7, verbose = FALSE)
  expect_equal(do.call(network_diffusion_with_pvalues, args),
               do.call(network_diffusion_with_pvalues, args))
})
