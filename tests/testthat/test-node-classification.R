# find_hubs() (high degree AND high betweenness) and find_bottlenecks() (low
# degree AND high betweenness) are near-mirror implementations that share a
# zscore/quantile thresholding scheme. Testing them together keeps the mirror
# honest.

g <- test_graph()

test_that("find_hubs() returns a complete node-level table", {
  res <- find_hubs(g, plot = FALSE)

  expect_named(res$result, c("node", "degree", "betweenness",
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
  res <- find_bottlenecks(g, plot = FALSE)

  expect_named(res$result, c("node", "degree", "betweenness",
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
  bn   <- find_bottlenecks(g, method = "quantile", plot = FALSE)

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
    hubs <- find_hubs(g, method = m, plot = FALSE)
    bn   <- find_bottlenecks(g, method = m, plot = FALSE)
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

test_that("an invalid thresholding method is rejected", {
  expect_error(find_hubs(g, method = "bogus", plot = FALSE), "should be one of")
  expect_error(find_bottlenecks(g, method = "bogus", plot = FALSE), "should be one of")
})
