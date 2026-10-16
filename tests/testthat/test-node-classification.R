# find_hubs() (high degree AND high betweenness) and find_bottlenecks() (low
# degree AND high betweenness) are near-mirror implementations that share a
# zscore/quantile thresholding scheme. Testing them together keeps the mirror
# honest.

g <- test_graph()

test_that("find_hubs() returns a complete node-level table", {
  res <- find_hubs(g, plot = FALSE)

  expect_named(res$result, c("node", "degree", "strength", "betweenness",
                             "degree_metric", "betweenness_metric", "is_hub"))
  expect_equal(nrow(res$result), igraph::vcount(g))
  expect_type(res$result$is_hub, "logical")
  expect_false(any(is.na(res$result$is_hub)))

  # Degree and betweenness are reported, not invented.
  expect_equal(res$result$degree,
               as.integer(igraph::degree(g)[res$result$node]),
               ignore_attr = TRUE)
})

test_that("find_bottlenecks() returns a complete node-level table", {
  res <- without_empty_selection(find_bottlenecks(g, plot = FALSE))

  expect_named(res$result, c("node", "degree", "strength", "betweenness",
                             "degree_metric", "betweenness_metric",
                             "is_bottleneck"))
  expect_equal(nrow(res$result), igraph::vcount(g))
  expect_type(res$result$is_bottleneck, "logical")
})

test_that("hubs sit at the top of the degree distribution", {
  res <- find_hubs(g, method = "quantile", degree_quantile = 0.9,
                   betweenness_quantile = 0.9, plot = FALSE)
  hubs <- res$result[res$result$is_hub, ]
  skip_if(nrow(hubs) == 0, "no hubs found in fixture")

  # A hub must out-degree the median node, by construction of the method.
  expect_gt(min(hubs$degree), stats::median(res$result$degree))
})

test_that("bottlenecks are low-degree but high-betweenness", {
  res <- find_bottlenecks(g, method = "quantile", degree_quantile = 0.75,
                          betweenness_quantile = 0.75, plot = FALSE)
  bn <- res$result[res$result$is_bottleneck, ]

  # Assert non-emptiness rather than skipping, so a fixture or igraph change
  # cannot turn the checks below into a vacuous pass.
  expect_gt(nrow(bn), 0)
  expect_lte(max(bn$degree), stats::quantile(res$result$degree, 0.75))
  expect_gte(min(bn$betweenness), stats::quantile(res$result$betweenness, 0.75))
})

hub_set <- function(res) res$result$node[res$result$is_hub]
bottleneck_set <- function(res) res$result$node[res$result$is_bottleneck]

test_that("hubs and bottlenecks are disjoint at their default thresholds", {
  # The defaults put the degree cutoffs far apart (0.95 for hubs, 0.25 for
  # bottlenecks), so the two sets cannot meet.
  hubs <- find_hubs(g, method = "quantile", plot = FALSE)
  bn   <- without_empty_selection(find_bottlenecks(g, method = "quantile",
                                                   plot = FALSE))

  expect_length(intersect(hub_set(hubs), bottleneck_set(bn)), 0L)
})

test_that("a shared degree cutoff puts boundary nodes in both sets", {
  # Both comparisons are inclusive -- hubs take degree >= cutoff, bottlenecks
  # degree <= cutoff -- so with one shared quantile every node sitting exactly on
  # the cutoff qualifies as both. Degrees are small integers, so that tie is
  # common rather than hypothetical. Documented so the inclusive boundary is a
  # deliberate choice rather than an accident.
  hubs <- find_hubs(g, method = "quantile", degree_quantile = 0.75,
                    betweenness_quantile = 0.75, plot = FALSE)
  bn   <- find_bottlenecks(g, method = "quantile", degree_quantile = 0.75,
                           betweenness_quantile = 0.75, plot = FALSE)

  both <- intersect(hub_set(hubs), bottleneck_set(bn))
  expect_gt(length(both), 0)

  # Every such node sits exactly on the shared degree cutoff.
  cutoff <- stats::quantile(hubs$result$degree_metric, 0.75, na.rm = TRUE)
  on_boundary <- hubs$result$degree_metric[match(both, hubs$result$node)]
  expect_true(all(on_boundary == cutoff))
})

test_that("both thresholding methods run and are reported in `method`", {
  for (m in c("zscore", "quantile")) {
    hubs <- without_empty_selection(find_hubs(g, method = m, plot = FALSE))
    bn   <- without_empty_selection(find_bottlenecks(g, method = m,
                                                     plot = FALSE))
    expect_match(hubs$method, m, fixed = TRUE)
    expect_match(bn$method, m, fixed = TRUE)
  }
})

test_that("stricter quantiles never enlarge the hub set", {
  loose  <- find_hubs(g, method = "quantile", degree_quantile = 0.80,
                      betweenness_quantile = 0.80, plot = FALSE)
  strict <- find_hubs(g, method = "quantile", degree_quantile = 0.95,
                      betweenness_quantile = 0.95, plot = FALSE)

  expect_lte(sum(strict$result$is_hub), sum(loose$result$is_hub))
  # Monotonicity: every strict hub is also a loose hub.
  expect_true(all(strict$result$node[strict$result$is_hub] %in%
                    loose$result$node[loose$result$is_hub]))
})

test_that("the is_hub vertex attribute matches the result table", {
  res <- find_hubs(g, plot = FALSE)
  from_attr  <- igraph::V(res$graph)$is_hub
  names(from_attr) <- igraph::V(res$graph)$name
  from_table <- stats::setNames(res$result$is_hub, res$result$node)

  expect_equal(from_attr[names(from_table)], from_table, ignore_attr = TRUE)
})

# --- Empty selections are reported, not returned in silence -----------------
#
# Both functions used to return a table with no selected rows and say nothing,
# which conflates three different situations: the graph has no such vertices,
# the threshold is unreachable, or the metric has no variance to threshold. The
# warning names which one occurred, because the remedy differs.

test_that("an unreachable z-threshold warns and says the threshold is the cause", {
  # A dense graph has a near-Gaussian degree distribution, so no vertex can
  # reach a z-score of 3 whatever the structure. This is the case reported from
  # the correlation-network analyses.
  set.seed(42)
  m <- matrix(stats::rnorm(80 * 40), nrow = 40)
  cm <- abs(stats::cor(m))
  diag(cm) <- 0
  dense <- igraph::graph_from_adjacency_matrix(
    cm > stats::quantile(cm[upper.tri(cm)], 0.7), mode = "undirected")
  igraph::V(dense)$name <- paste0("c", seq_len(igraph::vcount(dense)))

  # The property that makes this an unreachable threshold rather than an unmet
  # one: the best vertex in the graph still falls short.
  res <- without_empty_selection(find_hubs(dense, plot = FALSE))
  expect_false(any(res$result$is_hub))
  expect_lt(max(res$result$degree_metric), 3)

  expect_warning(find_hubs(dense, plot = FALSE), "found no hub nodes")
  expect_warning(find_hubs(dense, plot = FALSE), "no vertex reaches degree_threshold")
  # The message has to name the way out, not merely the problem.
  expect_warning(find_hubs(dense, plot = FALSE), "method = \"quantile\"",
                 fixed = TRUE)
})

test_that("an empty conjunction warns without blaming the thresholds", {
  # Both criteria are individually reachable here; no vertex meets both. That is
  # an answer about the graph, so the message must not suggest changing an
  # argument.
  #
  # method = "quantile" rather than the default: a quantile cutoff is by
  # construction a value some vertex attains, so "unreachable" is impossible and
  # the disjoint-conjunction branch is the only way the selection can be empty.
  # That makes this deterministic, where the previous version guarded itself with
  # two skip_if()s against a fixture that might not exercise the branch at all --
  # and after the default changed to "absolute" it skipped unconditionally,
  # because `degree_metric` is a raw neighbour count there and so never below -1.
  set.seed(7)
  dense <- igraph::sample_gnp(120, 0.25, directed = FALSE)
  igraph::V(dense)$name <- paste0("d", seq_len(igraph::vcount(dense)))

  res <- without_empty_selection(
    find_bottlenecks(dense, method = "quantile", plot = FALSE))
  expect_false(any(res$result$is_bottleneck))

  expect_warning(find_bottlenecks(dense, method = "quantile", plot = FALSE),
                 "no vertex passes both")
})

test_that("the absolute default also reports a disjoint conjunction", {
  # A star: the centre carries all the flow but has degree 4, and the leaves are
  # sparse but lie on no path. Each criterion is met, never by the same vertex.
  g <- igraph::make_star(5, mode = "undirected")
  igraph::V(g)$name <- paste0("s", seq_len(5))

  expect_warning(find_bottlenecks(g, plot = FALSE), "no vertex passes both")
})

test_that("a regular graph reports an undefined metric rather than NA rows", {
  # Every vertex of a ring has degree 2 and the same betweenness, so sd is 0 and
  # scale() returns NaN throughout. The classification used to come back NA for
  # every vertex and was then written onto the graph as an all-NA attribute,
  # which downstream filters silently dropped.
  ring <- ring_graph()

  expect_warning(find_hubs(ring, plot = FALSE), "same value at every vertex")
  expect_warning(find_bottlenecks(ring, plot = FALSE), "same value at every vertex")

  res <- without_empty_selection(find_hubs(ring, plot = FALSE))
  expect_false(any(is.na(res$result$is_hub)))
  expect_false(any(res$result$is_hub))
  # The z-score itself stays NA: it is genuinely undefined, and reporting 0
  # would claim every vertex sits exactly at the mean of a distribution that has
  # no spread.
  expect_true(all(is.na(res$result$degree_metric)))
  # The vertex attribute must agree with the table rather than carrying NA.
  expect_false(any(is.na(igraph::V(res$graph)$is_hub)))

  bn <- without_empty_selection(find_bottlenecks(ring, plot = FALSE))
  expect_false(any(is.na(bn$result$is_bottleneck)))
  expect_false(any(is.na(igraph::V(bn$graph)$is_bottleneck)))
})

test_that("a non-empty selection does not warn", {
  # The complement of the above: the warning must not fire merely because it
  # exists. Quantile thresholds on the scale-free fixture select something.
  expect_silent(find_hubs(g, method = "quantile", degree_quantile = 0.8,
                          betweenness_quantile = 0.8, plot = FALSE))
})

test_that("an invalid thresholding method is rejected", {
  expect_error(find_hubs(g, method = "bogus", plot = FALSE), "should be one of")
  expect_error(find_bottlenecks(g, method = "bogus", plot = FALSE), "should be one of")
})

# --- Zero betweenness is never high betweenness -----------------------------
#
# Both functions require high betweenness by definition. Betweenness is heavily
# zero-inflated on a sparse graph, so a quantile cutoff can land exactly on zero
# and `>=` then admits every vertex on that axis -- leaving a selection resting
# on degree alone while still claiming both criteria.

test_that("a vertex on no shortest path is not a bottleneck", {
  # A star is the exact case: the centre carries every path, all 20 leaves carry
  # none, so 95% of vertices have betweenness 0 and the 75th percentile is 0.
  # The leaves are also the lowest-degree vertices, so they passed the degree
  # criterion too -- the default call returned all 20 of them.
  star <- igraph::make_star(21, mode = "undirected", center = 1)
  igraph::V(star)$name <- paste0("s", seq_len(21))

  btw <- igraph::betweenness(star, normalized = TRUE)
  # The preconditions, so a future igraph change cannot make this vacuous.
  expect_equal(stats::quantile(log1p(btw), 0.75), 0, ignore_attr = TRUE)
  expect_equal(sum(btw == 0), 20L)

  res <- without_empty_selection(
    find_bottlenecks(star, method = "quantile", plot = FALSE))

  expect_false(any(res$result$is_bottleneck))
  # Specifically: the vertices that would have been selected are exactly the
  # ones with no betweenness at all.
  would_pass <- res$result$degree_metric <=
    stats::quantile(res$result$degree_metric, 0.25) &
    res$result$betweenness_metric >=
      stats::quantile(res$result$betweenness_metric, 0.75)
  expect_equal(sum(would_pass), 20L)
  expect_true(all(res$result$betweenness[would_pass] == 0))
})

test_that("a vertex on no shortest path is not a hub", {
  # The same guard on the mirror function. A complete graph has zero betweenness
  # everywhere, so no vertex can be a hub however well connected it is -- which
  # is the right answer: with every pair adjacent, no vertex mediates anything.
  full <- igraph::make_full_graph(12)
  igraph::V(full)$name <- paste0("f", seq_len(12))

  expect_true(all(igraph::betweenness(full) == 0))
  res <- without_empty_selection(
    find_hubs(full, method = "quantile", plot = FALSE))
  expect_false(any(res$result$is_hub))
})

test_that("a genuine bottleneck is still found", {
  # The complement: the guard must not remove the vertices the function exists
  # to find. The bridge vertices of a barbell carry every path between the two
  # halves, so their betweenness is far from zero.
  g2 <- igraph::disjoint_union(igraph::make_full_graph(6),
                               igraph::make_full_graph(6))
  igraph::V(g2)$name <- c(paste0("L", 1:6), paste0("R", 1:6))
  g2 <- igraph::add_vertices(g2, 2, name = c("b1", "b2"))
  g2 <- igraph::add_edges(g2, c("L1", "b1", "b1", "b2", "b2", "R1"))

  res <- find_bottlenecks(g2, plot = FALSE)
  expect_setequal(res$result$node[res$result$is_bottleneck], c("b1", "b2"))
  expect_true(all(res$result$betweenness[res$result$is_bottleneck] > 0))
})

test_that("the betweenness requirement does not change hubs on a scale-free graph", {
  # A high-degree vertex normally has paths through it, so the guard is a no-op
  # here. Asserted because "normally" was the previous contract and this makes
  # it an explicit one: the selection must come from the thresholds, not from
  # the guard.
  set.seed(91)
  pa <- igraph::sample_pa(300, power = 1.5, directed = FALSE)
  igraph::V(pa)$name <- paste0("v", seq_len(300))

  res <- without_empty_selection(
    find_hubs(pa, method = "quantile", plot = FALSE))
  expect_gt(sum(res$result$is_hub), 0)
  expect_true(all(res$result$betweenness[res$result$is_hub] > 0))
})

# --- find_bottlenecks(method = "absolute") ----------------------------------
#
# The default through 0.1.1 was "zscore", whose degree_threshold = -1 asked for
# low degree to be a statistical outlier. On a right-skewed degree distribution
# the left tail is short -- degree-1 vertices are the mode of a scale-free
# graph, not an anomaly -- so the documented default selected nothing at every
# size tested. "absolute" splits the two axes by the kind of claim each makes:
# low connectivity is absolute, high flow is relative.

test_that("the absolute default finds the bridge it is defined by", {
  # Two cliques joined by a two-vertex chain. b1 and b2 have degree 2 and carry
  # every path between the cliques: the low-degree, high-betweenness signature.
  g <- igraph::disjoint_union(igraph::make_full_graph(6), igraph::make_full_graph(6))
  igraph::V(g)$name <- c(paste0("L", 1:6), paste0("R", 1:6))
  g <- igraph::add_vertices(g, 2, name = c("b1", "b2"))
  g <- igraph::add_edges(g, c("L1", "b1", "b1", "b2", "b2", "R1"))

  res <- find_bottlenecks(g, plot = FALSE)
  expect_setequal(res$result$node[res$result$is_bottleneck], c("b1", "b2"))
  expect_match(res$method, "absolute")
})

test_that("the absolute default selects on a scale-free graph where zscore cannot", {
  # The point of the change. Not a claim that more is better -- a claim that the
  # old default was unreachable by construction, so it answered "no bottlenecks"
  # for every scale-free graph regardless of structure.
  set.seed(1)
  g <- igraph::sample_pa(2000, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(2000))

  res <- find_bottlenecks(g, plot = FALSE)
  expect_gt(sum(res$result$is_bottleneck), 0)
  expect_warning(find_bottlenecks(g, method = "zscore", plot = FALSE),
                 "no vertex reaches degree_threshold")
})

test_that("every absolutely-selected vertex is sparse and carries flow", {
  set.seed(2)
  g <- igraph::sample_pa(800, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(800))

  res <- find_bottlenecks(g, max_degree = 4, plot = FALSE)
  sel <- res$result[res$result$is_bottleneck, ]
  skip_if(nrow(sel) == 0, "fixture has no bottleneck at this cap")

  # The two criteria, asserted directly rather than via the count.
  expect_true(all(sel$degree <= 4))
  expect_true(all(sel$betweenness > 0))
})

test_that("max_degree is an absolute cap, so raising it cannot shrink the selection", {
  # Monotonicity: the degree axis is a cap on a fixed quantity, and the
  # betweenness cutoff does not depend on it, so the selection is nested.
  set.seed(3)
  g <- igraph::sample_pa(600, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(600))

  sel <- function(k) {
    r <- without_empty_selection(find_bottlenecks(g, max_degree = k, plot = FALSE))
    r$result$node[r$result$is_bottleneck]
  }
  s2 <- sel(2); s5 <- sel(5); s10 <- sel(10)
  expect_true(all(s2 %in% s5))
  expect_true(all(s5 %in% s10))
})

test_that("the degree cap is the neighbour count even under weights", {
  # An absolute threshold needs a unit comparable across graphs, and strength is
  # not one: scaling every weight by 1000 would otherwise empty the selection.
  g <- igraph::disjoint_union(igraph::make_full_graph(6), igraph::make_full_graph(6))
  igraph::V(g)$name <- c(paste0("L", 1:6), paste0("R", 1:6))
  g <- igraph::add_vertices(g, 2, name = c("b1", "b2"))
  g <- igraph::add_edges(g, c("L1", "b1", "b1", "b2", "b2", "R1"))

  igraph::E(g)$wt <- rep(2, igraph::ecount(g))
  small <- find_bottlenecks(g, weights = "wt", plot = FALSE)
  igraph::E(g)$wt <- rep(2000, igraph::ecount(g))
  large <- find_bottlenecks(g, weights = "wt", plot = FALSE)

  expect_setequal(small$result$node[small$result$is_bottleneck], c("b1", "b2"))
  expect_equal(small$result$node[small$result$is_bottleneck],
               large$result$node[large$result$is_bottleneck])
})

test_that("a regular graph selects nothing under the absolute default", {
  # The flaw this guard exists to prevent: "high betweenness" is a comparison,
  # so a quantile cutoff on a constant admits everything. Without the guard all
  # five vertices of a ring were returned as bottlenecks, which is the same
  # objection recorded against taking a fixed fraction of the graph.
  ring <- ring_graph()
  expect_warning(res <- find_bottlenecks(ring, plot = FALSE),
                 "same value at every vertex")
  expect_false(any(res$result$is_bottleneck))
  expect_false(any(is.na(res$result$is_bottleneck)))
  # Unlike the zscore path, nothing is standardized here, so the metric is
  # reported as measured rather than as NA.
  expect_false(any(is.na(res$result$betweenness_metric)))
})

test_that("a complete graph reports that no vertex carries any flow", {
  # Every pair is adjacent, so betweenness is zero throughout and the degree cap
  # is unreachable too. Both are named.
  res <- expect_warning(find_bottlenecks(igraph::make_full_graph(8), plot = FALSE),
                        "no vertex reaches")
  expect_false(any(res$result$is_bottleneck))
})

test_that("max_degree is validated eagerly", {
  g <- test_graph()
  expect_error(find_bottlenecks(g, max_degree = 0, plot = FALSE), "at least 1")
  expect_error(find_bottlenecks(g, max_degree = c(1, 2), plot = FALSE), "single finite")
  expect_error(find_bottlenecks(g, max_degree = NA, plot = FALSE), "single finite")
  expect_error(find_bottlenecks(g, max_degree = Inf, plot = FALSE), "single finite")
})

test_that("the older methods remain available and unchanged", {
  # The default changed; the alternatives did not.
  set.seed(4)
  g <- igraph::sample_pa(300, directed = FALSE)
  igraph::V(g)$name <- paste0("v", seq_len(300))

  q <- without_empty_selection(find_bottlenecks(g, method = "quantile", plot = FALSE))
  expect_match(q$method, "quantile")
  # The zero-betweenness guard still applies on every method.
  expect_true(all(q$result$betweenness[q$result$is_bottleneck] > 0))

  z <- without_empty_selection(find_bottlenecks(g, method = "zscore", plot = FALSE))
  expect_match(z$method, "zscore")
})
