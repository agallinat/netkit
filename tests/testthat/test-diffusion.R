# The diffusion cluster is netkit's most intricate code, and the
# adjacency-normalise -> Laplacian -> Cholesky/transition-matrix block is
# currently implemented twice: once inside network_diffusion() and once inside
# prepare_diffusion(). The precompute-equivalence test below is the guard that
# makes consolidating those two copies safe.

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

test_that("precompute() reproduces the un-precomputed result exactly", {
  # If this ever fails, the two copies of the kernel-construction block have
  # drifted apart.
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
