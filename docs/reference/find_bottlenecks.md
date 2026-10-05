# Identify and Bottleneck Nodes in a Network

Identifies bottleneck nodes in an `igraph` network as those with low
degree and high betweenness centrality. The function supports both
standardized (z-score) and quantile-based thresholding. Optionally, it
produces a 2D scatter plot with bottlenecks highlighted.

## Usage

``` r
find_bottlenecks(
  graph,
  method = c("zscore", "quantile"),
  degree_threshold = -1,
  betweenness_threshold = 1,
  degree_quantile = 0.25,
  betweenness_quantile = 0.75,
  log_transform = TRUE,
  plot = TRUE,
  focus_color = "skyblue",
  bottleneck_names = TRUE,
  bottleneck_cex = 3,
  gg_extra = list(),
  weights = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph:

  An `igraph` object representing the network to analyze or a data frame
  containing a symbolic edge list in the first two columns. Additional
  columns are considered as edge attributes.

- method:

  Character. Method to define bottlenecks: `"zscore"` or `"quantile"`.

- degree_threshold:

  Numeric. Upper threshold for standardized degree (only used if
  `method = "zscore"`).

- betweenness_threshold:

  Numeric. Lower threshold for standardized betweenness (used in both
  methods).

- degree_quantile:

  Numeric between 0 and 1. Quantile threshold for degree (used if
  `method = "quantile"`).

- betweenness_quantile:

  Numeric between 0 and 1. Quantile threshold for betweenness (used if
  `method = "quantile"`).

- log_transform:

  Logical. If `TRUE`, applies `log1p` transformation to degree and
  betweenness.

- plot:

  Logical. If `TRUE`, generates a plot of degree vs. betweenness
  highlighting bottlenecks.

- focus_color:

  Character. Color to display in the focus area of the plot (bottlenecks
  region).

- bottleneck_names:

  Logical. If `TRUE`, labels bottleneck nodes on the plot.

- bottleneck_cex:

  Numeric. Font size scaling for bottleneck labels on the plot.

- gg_extra:

  List. Additional user-defined layers for the returned ggplot. eg.
  list(ylim(-2,2), theme_bw(), theme(legend.position = "none"))

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  When supplied, degree becomes vertex *strength* and betweenness is
  computed on edge *costs*. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

## Value

A list with the following components:

- `plot`:

  A scatter plot of degree vs. betweenness highlighting bottlenecks, or
  `NULL` when `plot = FALSE`. The element is always present, so the
  return shape does not depend on the arguments. Note this is a
  `ggExtraPlot` (a grid gtable), not a `ggplot`: the scatterplot is
  passed through
  [`ggMarginal()`](https://rdrr.io/pkg/ggExtra/man/ggMarginal.html) to
  add the marginal histograms, which assembles it, so it can be printed
  but not extended with `+`. To customize it, pass layers through
  `gg_extra`; those are added before the assembly step.

- `method`:

  A message describing the method and thresholds used.

- `result`:

  A `tibble` with node name, degree, betweenness, transformed metrics,
  and bottleneck status.

- `graph`:

  The original graph with a new vertex attribute `is_bottleneck`
  (logical).

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
g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
res <- find_bottlenecks(g, method = "quantile", plot = FALSE)
res$method
#> [1] "Bottlenecks identified by method: quantile with Degree metric threshold = 0.693147180559945 and Betweenness metric threshold = 0.0250013022054173 (unweighted)"
head(res$result)
#> # A tibble: 6 × 7
#>   node  degree strength betweenness degree_metric betweenness_metric
#>   <chr>  <dbl>    <dbl>       <dbl>         <dbl>              <dbl>
#> 1 1         28       28      0.951          3.37              0.668 
#> 2 2          5        5      0.123          1.79              0.116 
#> 3 3          3        3      0.0987         1.39              0.0941
#> 4 4          4        4      0.0750         1.61              0.0723
#> 5 5          2        2      0.0253         1.10              0.0250
#> 6 6          1        1      0              0.693             0     
#> # ℹ 1 more variable: is_bottleneck <lgl>
```
