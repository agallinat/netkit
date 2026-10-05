# Plot Complementary Cumulative Degree Distribution (CCDF)

This function plots the complementary cumulative distribution function
(CCDF) of node degrees in a network and optionally overlays power-law
reference curves.

## Usage

``` r
plot_CCDF(
  graph,
  keep_direction = TRUE,
  remove_singles = FALSE,
  show_PL = TRUE,
  PL_exponents = c(2, 3),
  colors = c("#000831", "#e41a1c", "darkgreen", "#9c52f2", "#b8b8ff"),
  label.size = 12,
  weights = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph:

  An `igraph` object representing the network to analyze or a data frame
  containing a symbolic edge list in the first two columns. Additional
  columns are considered as edge attributes.

- keep_direction:

  Logical. Only for directed graphs. If `TRUE`, CCDF curves are drawn
  for 'in'-degree, 'out'-degree, and 'all'-degree distributions. `FALSE`
  to ignore directionality.

- remove_singles:

  Logical. If `TRUE`, nodes with degree 0 are removed from the graph
  before computing the CCDF. Default is `FALSE`.

- show_PL:

  Logical. If `TRUE`, overlays theoretical power-law reference lines of
  the form \\P(K \> k) \sim k^{-\gamma}\\. Default is `TRUE`.

- PL_exponents:

  Numeric vector. The \\\gamma\\ exponents for the power-law curves.
  Default is `c(2, 3)`.

- colors:

  Optional character vector. Custom colors for the graph curve and
  power-law lines. If `NULL`, default colors are used.

- label.size:

  Numeric. Font size for axis labels and theme. Passed to
  `theme_minimal`.

- weights:

  Optional edge weights: `NULL` (default) to plot the degree
  distribution, or the name of an edge attribute / a numeric vector of
  length `igraph::ecount(graph)` to plot the vertex *strength*
  distribution instead. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

## Value

A `ggplot2` object showing the CCDF of node degrees on a log-log scale.

## Edge weights

Functions that can use edge weights take two arguments:

- `weights`:

  `NULL` (the default) to ignore edge weights; the name of an edge
  attribute, such as `"weight"`; or a numeric vector with one value per
  edge, in `igraph::E(graph)` order.

- `weight_type`:

  `"strength"` (the default) if a larger value means a more tightly
  connected pair – confidence scores, correlations, co-expression, read
  counts, interaction scores. `"distance"` if a larger value means
  further apart – costs, dissimilarities, reaction times.

Declaring which you have is not bookkeeping. igraph reads the `weight`
attribute implicitly and gives it *opposite* meanings in different
functions: a cost in
[`igraph::betweenness()`](https://r.igraph.org/reference/betweenness.html),
[`igraph::distances()`](https://r.igraph.org/reference/distances.html),
[`igraph::diameter()`](https://r.igraph.org/reference/diameter.html) and
[`igraph::mean_distance()`](https://r.igraph.org/reference/distances.html),
but a strength in
[`igraph::cluster_louvain()`](https://r.igraph.org/reference/cluster_louvain.html)
and the other community detection algorithms. Attaching a confidence
score and letting that happen implicitly therefore inverts every
path-based metric – a high-confidence interaction is treated as a long
distance – while community detection reads the same numbers the way you
intended.

netkit resolves `weights` and `weight_type` once per call and derives
both a strength and a distance vector from them, so each metric receives
the one it needs. A `"strength"` is converted to a distance by
reciprocal (\\1/w\\); a `"distance"` is converted to a strength by
reflection (\\\max(w) - w + \min(w)\\), which keeps a zero distance
finite.

## Examples

``` r
g <- igraph::sample_pa(200, power = 1.5, directed = FALSE)
plot_CCDF(g)


# Compare against reference power-law slopes.
plot_CCDF(g, PL_exponents = c(2, 2.5, 3))

```
