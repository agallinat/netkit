# calculate_roles() implements the Guimera & Amaral (2005) R1-R7 cartography.
# Until recently every call errored with `could not find function "neighbors"`,
# because igraph::neighbors was called unqualified and was absent from
# NAMESPACE. The first test is the regression guard for that.

g <- test_graph()

test_that("calculate_roles() runs without igraph attached", {
  # The original bug only manifested when igraph was *not* on the search path,
  # which is exactly how the package is used via netkit::. Guard that
  # precondition explicitly so the test cannot pass for the wrong reason.
  skip_if("package:igraph" %in% search(),
          "igraph is attached; this test needs it off the search path")

  expect_no_error(calculate_roles(g, cluster.method = "louvain", plot = FALSE))
})

test_that("the result table has exactly the documented columns", {
  res <- calculate_roles(g, cluster.method = "louvain", plot = FALSE)

  # A stray `stringsAsFactors` column used to leak in here: tibble() has no such
  # argument, so it was silently recycled into a data column.
  expect_named(res$result, c("node", "module", "z", "p", "role"))
  expect_type(res$result$node, "character")
  expect_type(res$result$module, "integer")
  expect_type(res$result$z, "double")
  expect_type(res$result$p, "double")
  expect_type(res$result$role, "character")
})

test_that("roles_definitions describes all seven roles", {
  res <- calculate_roles(g, cluster.method = "louvain", plot = FALSE)

  expect_s3_class(res$roles_definitions, "data.frame")
  expect_equal(nrow(res$roles_definitions), 7L)
  expect_named(res$roles_definitions, c("Name", "Description", "Condition"))
  expect_equal(res$roles_definitions$Name, paste0("R", 1:7))
})

test_that("z and P fall in their valid ranges", {
  res <- calculate_roles(g, cluster.method = "louvain", plot = FALSE)

  # The participation coefficient is bounded on [0, 1] by definition.
  expect_true(all(res$result$p >= 0 & res$result$p <= 1))
  expect_false(any(is.na(res$result$p)))
  expect_false(any(is.na(res$result$z)))
})

test_that("assigned roles are consistent with the z/P thresholds", {
  res <- calculate_roles(g, cluster.method = "louvain", plot = FALSE, hub_z = 2.5)
  r <- res$result

  expect_true(all(r$role %in% paste0("R", 1:7)))

  # Non-hubs (z < hub_z) get R1-R4; hubs get R5-R7.
  expect_true(all(r$role[r$z <  2.5] %in% c("R1", "R2", "R3", "R4")))
  expect_true(all(r$role[r$z >= 2.5] %in% c("R5", "R6", "R7")))

  # Spot-check the P cut points on the non-hub branch.
  expect_true(all(r$p[r$role == "R1"] <= 0.05))
  expect_true(all(r$p[r$role == "R2"] >  0.05 & r$p[r$role == "R2"] <= 0.60))
  expect_true(all(r$p[r$role == "R3"] >  0.60 & r$p[r$role == "R3"] <= 0.80))
  expect_true(all(r$p[r$role == "R4"] >  0.80))
})

test_that("hub_z shifts the hub/non-hub boundary", {
  strict <- calculate_roles(g, cluster.method = "louvain", plot = FALSE, hub_z = 2.5)
  loose  <- calculate_roles(g, cluster.method = "louvain", plot = FALSE, hub_z = 0.5)

  n_hubs <- function(res) sum(res$result$role %in% c("R5", "R6", "R7"))
  expect_gte(n_hubs(loose), n_hubs(strict))

  # The thresholds are echoed back in roles_definitions.
  expect_match(loose$roles_definitions$Condition[1], "0.5", fixed = TRUE)
})

test_that("a precomputed membership vector is accepted", {
  set.seed(1)
  comm <- igraph::cluster_louvain(g)

  from_object <- calculate_roles(g, communities = comm, plot = FALSE)
  expect_equal(nrow(from_object$result), igraph::vcount(g))

  from_vector <- calculate_roles(g, communities = igraph::membership(comm),
                                 plot = FALSE)
  expect_equal(nrow(from_vector$result), igraph::vcount(g))
  expect_equal(from_object$result$role, from_vector$result$role)
})

test_that("a malformed communities argument is rejected", {
  expect_error(calculate_roles(g, communities = 1:3, plot = FALSE),
               "clustering object, a membership vector or NULL")
})

test_that("falling back to find_modules() can drop nodes from small modules", {
  # Documented quirk rather than endorsed behaviour: with communities = NULL the
  # function delegates to find_modules(), whose default min_size = 3 discards
  # small modules, so those nodes never reach the roles table.
  res <- calculate_roles(g, cluster.method = "louvain", plot = FALSE)

  expect_lte(nrow(res$result), igraph::vcount(g))
  expect_true(all(res$result$node %in% igraph::V(g)$name))
  expect_false(any(duplicated(res$result$node)))
})
