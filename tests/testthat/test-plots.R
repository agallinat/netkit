# netkit builds diagnostics as ggplot objects and returns them rather than
# printing, with two deliberate exceptions: plot_Net() and highlight_nodes(),
# which draw with base plot.igraph because they manipulate par() for the colour
# legend. These tests pin that split.

g <- test_graph()

test_that("plot_CCDF() returns a ggplot without drawing", {
  p <- plot_CCDF(g)
  expect_s3_class(p, "ggplot")
  # A returned-but-not-printed plot must survive being built.
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_CCDF() handles directed graphs and the keep_direction switch", {
  dg <- test_graph_directed()
  expect_s3_class(plot_CCDF(dg, keep_direction = TRUE), "ggplot")
  expect_s3_class(plot_CCDF(dg, keep_direction = FALSE), "ggplot")
})

test_that("plot_CCDF() options are accepted", {
  expect_s3_class(plot_CCDF(g, show_PL = FALSE), "ggplot")
  expect_s3_class(plot_CCDF(g, remove_singles = TRUE), "ggplot")
  expect_s3_class(plot_CCDF(g, PL_exponents = c(2, 2.5, 3)), "ggplot")
})

test_that("plot_Net() draws and returns nothing", {
  expect_null(draw_quietly(plot_Net(g)))
  expect_null(draw_quietly(plot_Net(g, label = TRUE)))
  expect_null(draw_quietly(plot_Net(g, node.degree.map = FALSE)))
  # Mapping a numeric vertex attribute to colour.
  expect_null(draw_quietly(plot_Net(g, color = "score", node.degree.map = FALSE)))
})

test_that("plot_Net() accepts a supplied layout", {
  lay <- layout_horizontal_tree(g)
  expect_null(draw_quietly(plot_Net(g, layout = lay)))
})

test_that("highlight_nodes() draws for every method, and combinations", {
  nodes <- igraph::V(g)$name[1:5]
  for (m in c("label", "fill", "outline")) {
    expect_null(draw_quietly(highlight_nodes(g, nodes, method = m)))
  }
  # method is a set, not a single choice: the default is all three.
  expect_null(draw_quietly(highlight_nodes(g, nodes)))
  expect_null(draw_quietly(highlight_nodes(g, nodes, method = c("label", "fill"))))
})

test_that("highlight_nodes() does not warn about NA frame colours", {
  # Methods other than "outline" used to set frame.color to NA for every vertex,
  # which made plot.igraph warn "vertex attribute frame.color contains NAs.
  # Replacing with default value black". plot_Net() already defaults the frame to
  # the fill colour when none is given, so the attribute is simply left unset.
  nodes <- igraph::V(g)$name[1:5]
  for (m in list("label", "fill", "outline", c("label", "fill"))) {
    expect_no_warning(draw_quietly(highlight_nodes(g, nodes, method = m)))
  }
})

test_that("a stale frame.color attribute does not leak into a non-outline call", {
  nodes <- igraph::V(g)$name[1:5]
  gf <- g
  igraph::V(gf)$frame.color <- NA   # as a previous outline-less call used to leave it
  expect_no_warning(draw_quietly(highlight_nodes(gf, nodes, method = "fill")))
})

test_that("layout_horizontal_tree() returns one coordinate pair per vertex", {
  lay <- layout_horizontal_tree(g)
  expect_true(is.matrix(lay))
  expect_equal(dim(lay), c(igraph::vcount(g), 2L))
  expect_true(all(is.finite(lay)))
})

test_that("layout_horizontal_tree() rotates igraph's vertical tree layout", {
  # The function post-rotates layout_as_tree(), so the x/y spans should swap
  # relative to the unrotated layout.
  vertical <- igraph::layout_as_tree(g)
  horizontal <- layout_horizontal_tree(g)

  span <- function(m, i) diff(range(m[, i]))
  expect_equal(span(horizontal, 1), span(vertical, 2), tolerance = 1e-8)
  expect_equal(span(horizontal, 2), span(vertical, 1), tolerance = 1e-8)
})

test_that("frame.color is resolved from frame.color, not from color", {
  # The branch deciding whether frame.color names a vertex attribute used to
  # test `color` instead, so the two arguments could only be used in matching
  # forms. Each case below took the wrong branch: the first looked up an
  # attribute named "black", the second painted every frame with the literal
  # string "grp". Both are asserted because they are different wrong branches.
  g <- test_graph(n = 20)
  igraph::V(g)$grp <- rep(c("red", "blue"), length.out = igraph::vcount(g))

  # Attribute-mapped fill with a literal frame.
  expect_no_error(draw_quietly(plot_Net(g, color = "grp", frame.color = "black")))

  # Literal fill with an attribute-mapped frame. Under the old branch this
  # passed the literal string "grp" to plot.igraph as a colour and errored.
  expect_no_error(draw_quietly(plot_Net(g, color = "skyblue",
                                        frame.color = "grp")))

  # plot_Net() draws rather than returning the graph, so the resolution cannot
  # be inspected directly. This discriminates instead: a literal frame.color
  # that is not a colour name must reach plot.igraph and fail there. Under the
  # old branch, `color` naming an attribute sent this down the attribute path,
  # which looked up a non-existent attribute, got NULL, set nothing, and drew
  # happily -- so the argument was silently discarded.
  expect_error(draw_quietly(plot_Net(g, color = "grp",
                                     frame.color = "notacolour")))
})
