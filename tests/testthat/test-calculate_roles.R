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

# --- z and P against a direct reference implementation ----------------------
#
# z and P are computed by whole-graph aggregation: within-module degree from one
# pass over the edge list, and the per-module neighbor tally from one sparse
# matrix. The implementations they replaced built an induced_subgraph() per
# module and called neighbors() once per vertex. These tests pin the new
# arithmetic against that older, slower, obviously-correct formulation rather
# than against recorded output, so a future rewrite has something to answer to.

ref_z <- function(g, memb, wname = NULL) {
  gn <- igraph::V(g)$name
  if (is.null(names(memb))) memb <- stats::setNames(memb, gn)
  z <- stats::setNames(rep(NA_real_, length(gn)), gn)
  for (mod in unique(memb[!is.na(memb)])) {
    nodes <- names(memb)[!is.na(memb) & memb == mod]
    if (length(nodes) <= 1) {
      z[nodes] <- 0
      next
    }
    gg <- if (is.null(wname)) {
      g
    } else {
      igraph::set_edge_attr(g, ".s", value = igraph::edge_attr(g, wname))
    }
    sub <- igraph::induced_subgraph(gg, nodes)
    ki <- if (is.null(wname)) {
      igraph::degree(sub)
    } else {
      igraph::strength(sub, weights = igraph::edge_attr(sub, ".s"))
    }
    s <- stats::sd(ki)
    z[names(ki)] <- if (is.na(s) || s == 0) rep(0, length(ki)) else (ki - mean(ki)) / s
  }
  z
}

ref_p <- function(g, memb, wname = NULL) {
  gn <- igraph::V(g)$name
  if (is.null(names(memb))) memb <- stats::setNames(memb, gn)
  strength <- if (is.null(wname)) NULL else igraph::edge_attr(g, wname)
  vapply(gn, function(node) {
    nbrs <- igraph::neighbors(g, node, mode = "all")
    if (length(nbrs) == 0) return(0)
    if (is.null(strength)) {
      mods <- memb[gn[as.integer(nbrs)]]
      mods <- mods[!is.na(mods)]
      if (!length(mods)) return(NA_real_)
      k <- as.numeric(table(mods))
    } else {
      e <- igraph::incident(g, node, mode = "all")
      en <- igraph::ends(g, e, names = TRUE)
      far <- ifelse(en[, 1] == node, en[, 2], en[, 1])
      m <- memb[far]
      ok <- !is.na(m)
      if (!any(ok)) return(NA_real_)
      k <- as.numeric(tapply(strength[as.integer(e)][ok], m[ok], sum))
    }
    1 - sum((k / sum(k))^2)
  }, numeric(1))
}

expect_matches_reference <- function(g, memb, wname = NULL) {
  res <- suppressWarnings(
    calculate_roles(g, communities = memb, plot = FALSE, weights = wname))
  got_z <- stats::setNames(res$result$z, res$result$node)
  got_p <- stats::setNames(res$result$p, res$result$node)
  expect_equal(unname(got_z), unname(ref_z(g, memb, wname)[names(got_z)]))
  expect_equal(unname(got_p), unname(ref_p(g, memb, wname)[names(got_p)]))
}

test_that("z and P match a per-module, per-vertex reference implementation", {
  set.seed(9)
  g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
  igraph::V(g)$name <- paste0("n", seq_len(60))
  expect_matches_reference(g, igraph::membership(igraph::cluster_louvain(g)))

  igraph::E(g)$wt <- stats::runif(igraph::ecount(g), 0.5, 5)
  expect_matches_reference(g, igraph::membership(igraph::cluster_louvain(g)), "wt")
})

test_that("z and P match the reference on a directed graph", {
  # Within-module degree and the neighbor tally both read mode = "all", so an
  # edge contributes to each endpoint regardless of direction.
  set.seed(4)
  g <- igraph::sample_pa(50, directed = TRUE)
  igraph::V(g)$name <- paste0("n", seq_len(50))
  memb <- stats::setNames(rep(1:5, length.out = 50), igraph::V(g)$name)

  expect_matches_reference(g, memb)
  igraph::E(g)$wt <- stats::runif(igraph::ecount(g), 0.5, 5)
  expect_matches_reference(g, memb, "wt")
})

test_that("self-loops and parallel edges are tallied as igraph counts them", {
  # The edge-list pass has to reproduce degree(loops = TRUE): a self-loop adds
  # its contribution twice, once per endpoint, and parallel edges count
  # separately. Getting this wrong would shift z without making it invalid.
  g <- igraph::graph_from_edgelist(
    as.matrix(data.frame(a = c(1, 1, 2, 2, 3, 4, 4),
                         b = c(2, 2, 3, 2, 1, 4, 5))), directed = FALSE)
  igraph::V(g)$name <- paste0("n", seq_len(igraph::vcount(g)))
  memb <- stats::setNames(c(1, 1, 2, 2, 2), igraph::V(g)$name)

  expect_matches_reference(g, memb)
  igraph::E(g)$wt <- c(1, 2, 3, 4, 5, 6, 7)
  expect_matches_reference(g, memb, "wt")
})

test_that("an edge to an unknown-module neighbour does not corrupt the tally", {
  # `TRUE & NA` is NA. Guarding only one endpoint when testing whether an edge is
  # internal to a module put NA into a subscript, which errors outright -- but the
  # same slip in the P tally would have silently mixed an unknown module into the
  # denominator. Both ends are guarded, and the NA-module vertices keep NA roles
  # while their neighbours stay computable.
  set.seed(11)
  g <- igraph::add_vertices(igraph::sample_gnp(40, 0.1, directed = FALSE), 2)
  igraph::V(g)$name <- paste0("n", seq_len(igraph::vcount(g)))
  memb <- as.integer(igraph::membership(igraph::cluster_louvain(g)))
  names(memb) <- igraph::V(g)$name
  memb[c(3, 7, 11, 19, 23, 31, 37, 41)] <- NA

  expect_matches_reference(g, memb)
  igraph::E(g)$wt <- stats::runif(igraph::ecount(g), 0.5, 5)
  expect_matches_reference(g, memb, "wt")

  res <- suppressWarnings(calculate_roles(g, communities = memb, plot = FALSE))
  # The unknown-module vertices are reported, not quietly given a module.
  expect_true(all(is.na(res$result$role[is.na(res$result$module)])))
  # ... and vertices that do have a module still get a role, so one NA neighbour
  # does not propagate across the graph.
  known <- !is.na(res$result$module) & res$result$module %in% res$result$module
  expect_true(any(!is.na(res$result$role[known])))
})

test_that("every vertex in a singleton module gets z = 0", {
  # A module of one has no within-module variance to standardize against. The
  # per-module sd is NA there rather than 0, so the zero-variance branch has to
  # catch both.
  set.seed(5)
  g <- igraph::sample_gnp(80, 0.06, directed = FALSE)
  igraph::V(g)$name <- paste0("n", seq_len(80))
  memb <- stats::setNames(seq_len(80), igraph::V(g)$name)

  res <- suppressWarnings(calculate_roles(g, communities = memb, plot = FALSE))
  expect_true(all(res$result$z == 0))
  expect_matches_reference(g, memb)
})

test_that("an edgeless graph yields z = 0 and P = 0 throughout", {
  g <- igraph::make_empty_graph(5, directed = FALSE)
  igraph::V(g)$name <- paste0("n", seq_len(5))
  res <- calculate_roles(g, communities = stats::setNames(rep(1L, 5), igraph::V(g)$name),
                         plot = FALSE)
  expect_true(all(res$result$z == 0))
  # No neighbours at all is P = 0, distinct from neighbours of unknown module,
  # which is P = NA.
  expect_true(all(res$result$p == 0))
})

test_that("the result table carries plain unnamed numeric columns", {
  # lengths() and rowsum() both carry names through, and ifelse() propagates the
  # attributes of its first argument, so the aggregation can hand back a named
  # vector that would surface as named tibble columns.
  g <- test_graph()
  res <- suppressWarnings(calculate_roles(g, cluster.method = "louvain", plot = FALSE))
  expect_null(names(res$result$z))
  expect_null(names(res$result$p))
  expect_null(names(res$result$module))
})
