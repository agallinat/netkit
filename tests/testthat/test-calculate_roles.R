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

test_that("every vertex receives a role", {
  # find_modules() is called with min_size = 1 so that no module is discarded.
  # Previously it ran at the default min_size = 3 and nodes in small modules never
  # reached the roles table.
  for (m in c("louvain", "walktrap", "infomap")) {
    set.seed(7)
    res <- calculate_roles(g, cluster.method = m, plot = FALSE)
    expect_setequal(res$result$node, igraph::V(g)$name)
    expect_false(any(duplicated(res$result$node)))
  }
})

test_that("the participation coefficient is unaffected by module-size filtering", {
  # Regression guard for a silent numerical bug. When a module was discarded, its
  # nodes' membership became NA; table() dropped those NAs from a neighbor's tally
  # while the node's full degree was still used as the denominator, so P came out
  # too high for *retained* nodes. On a 150-node graph with walktrap this gave 18
  # wrong coefficients and 7 wrong roles.
  set.seed(3)
  gw <- igraph::sample_gnp(150, 0.02, directed = FALSE)
  igraph::V(gw)$name <- paste0("b", seq_len(150))

  set.seed(7)
  via_find_modules <- calculate_roles(gw, cluster.method = "walktrap", plot = FALSE)
  set.seed(7)
  via_membership <- calculate_roles(
    gw, communities = igraph::membership(igraph::cluster_walktrap(gw)), plot = FALSE
  )

  expect_equal(nrow(via_find_modules$result), igraph::vcount(gw))

  a <- via_find_modules$result[order(via_find_modules$result$node), ]
  b <- via_membership$result[order(via_membership$result$node), ]
  expect_equal(a$p, b$p)
  expect_equal(a$role, b$role)
})

test_that("the participation coefficient stays within bounds for pendant nodes", {
  # A degree-1 node has all its edges inside one module, so P must be exactly 0.
  star <- igraph::make_star(10, mode = "undirected")
  igraph::V(star)$name <- paste0("s", seq_len(10))
  res <- calculate_roles(star, communities = rep(1L, 10), plot = FALSE)

  leaves <- res$result[res$result$node != "s1", ]
  expect_true(all(leaves$p == 0))
  expect_true(all(res$result$p >= 0 & res$result$p <= 1))
})

test_that("vertices with no module assignment are reported, not dropped silently", {
  # spinglass cannot run on a disconnected graph, so find_modules() restricts it to
  # the largest connected component. The resulting partial coverage now warns.
  set.seed(42)
  gd <- igraph::sample_gnp(80, 0.02, directed = FALSE)
  igraph::V(gd)$name <- paste0("s", seq_len(80))
  expect_gt(igraph::components(gd)$no, 1)

  # Two warnings are expected and both are asserted: find_modules() reports the
  # fallback to the LCC, and calculate_roles() reports the resulting partial
  # coverage.
  expect_warning(
    expect_warning(
      res <- calculate_roles(gd, cluster.method = "spinglass", plot = FALSE),
      "vertices have no module assignment"
    ),
    "cannot work with unconnected graph"
  )
  expect_lt(nrow(res$result), igraph::vcount(gd))
})

test_that("a membership vector containing NAs is handled rather than erroring", {
  # Previously this failed with "Invalid vertex names": the NA module was treated as
  # a module, and induced_subgraph() was handed NA node names.
  set.seed(1)
  gn <- igraph::sample_pa(40, directed = FALSE)
  igraph::V(gn)$name <- paste0("n", seq_len(40))
  m <- igraph::membership(igraph::cluster_louvain(gn))
  m[1:5] <- NA

  res <- suppressWarnings(calculate_roles(gn, communities = m, plot = FALSE))

  expect_equal(nrow(res$result), 40L)
  expect_equal(sum(is.na(res$result$module)), 5L)
  # Undefined inputs propagate to an undefined role rather than a wrong one.
  expect_true(all(is.na(res$result$role[is.na(res$result$module)])))
  expect_true(all(res$result$p >= 0 & res$result$p <= 1, na.rm = TRUE))
})

# --- Regression: one source of truth for the role thresholds ----------------
#
# The five participation-coefficient boundaries used to be written out three
# times -- in the classifier, in the roles_definitions table handed to the
# caller, and in the shaded bands of the plot -- and had drifted: the R5/R6
# boundary was 0.30 in the classifier but 0.25 in the other two. So the shaded
# "provincial hub" region disagreed with the classification it illustrated, and
# the table documented a threshold the code did not use. Both now read the same
# vector, and these tests pin that they agree.

test_that("roles_definitions reports the boundaries the classifier actually used", {
  res <- calculate_roles(test_graph(), cluster.method = "louvain", plot = FALSE)
  conds <- res$roles_definitions$Condition

  # The published hub boundary is 0.30, and this is the one that had drifted.
  expect_match(conds[5], "P <= 0.3", fixed = TRUE)
  expect_match(conds[6], "0.3 < P", fixed = TRUE)
  expect_false(any(grepl("0.25", conds, fixed = TRUE)))

  # Non-hub boundaries, as published in Guimera & Amaral (2005).
  expect_match(conds[1], "P <= 0.05", fixed = TRUE)
  expect_match(conds[2], "0.05 < P & P <= 0.62", fixed = TRUE)
  expect_match(conds[3], "0.62 < P & P <= 0.8", fixed = TRUE)
})

test_that("the classifier honours overridden thresholds", {
  g <- test_graph()

  # Push every non-hub node into R1 by moving the R1/R2 boundary near its
  # maximum. The boundaries must stay strictly increasing for each role to
  # remain reachable, hence 0.98/0.99/1 rather than three copies of 1.
  res <- calculate_roles(g, cluster.method = "louvain", plot = FALSE,
                         thresholds = c(R1_R2 = 0.98, R2_R3 = 0.99, R3_R4 = 1))
  non_hub <- res$result[!is.na(res$result$z) & res$result$z < 2.5 &
                          !is.na(res$result$p), ]
  expect_true(all(non_hub$role == "R1"))

  # And that the reported definitions move with them.
  expect_match(res$roles_definitions$Condition[1], "P <= 0.98", fixed = TRUE)
})

test_that("invalid thresholds are rejected eagerly", {
  g <- test_graph()
  expect_error(calculate_roles(g, cluster.method = "louvain", plot = FALSE,
                               thresholds = c(nonsense = 0.5)),
               "Unknown 'thresholds' name")
  expect_error(calculate_roles(g, cluster.method = "louvain", plot = FALSE,
                               thresholds = c(R1_R2 = 1.5)),
               "within \\[0, 1\\]")
  expect_error(calculate_roles(g, cluster.method = "louvain", plot = FALSE,
                               thresholds = c(R1_R2 = 0.9, R2_R3 = 0.1)),
               "unreachable")
  expect_error(calculate_roles(g, cluster.method = "louvain", plot = FALSE,
                               thresholds = 0.5),
               "named numeric vector")
})

# --- Regression: the return contract ----------------------------------------

test_that("calculate_roles() returns the full plot/result/graph/method shape", {
  res <- calculate_roles(test_graph(), cluster.method = "louvain", plot = FALSE)
  expect_true(all(c("plot", "result", "graph", "method") %in% names(res)))
  expect_s3_class(res$graph, "igraph")
  expect_type(res$method, "character")
  expect_length(res$method, 1)
})

test_that("the returned graph carries the classification and chains onward", {
  g <- test_graph()
  res <- calculate_roles(g, cluster.method = "louvain", plot = FALSE)

  expect_true(all(c("module", "role_z", "role_p", "role") %in%
                    igraph::vertex_attr_names(res$graph)))
  expect_type(igraph::V(res$graph)$role, "character")

  # Attributes are keyed by name, not by position: the graph attribute must
  # agree with the result table row for the same node, not merely be the same
  # length. Shuffling the table would break a positional assignment.
  ord <- match(igraph::V(res$graph)$name, res$result$node)
  expect_equal(igraph::V(res$graph)$role, res$result$role[ord])
  expect_equal(igraph::V(res$graph)$role_p, res$result$p[ord])

  # Pre-existing attributes survive, so the graph can keep chaining.
  expect_true("score" %in% igraph::vertex_attr_names(res$graph))
  expect_silent(find_hubs(res$graph, plot = FALSE))
})
