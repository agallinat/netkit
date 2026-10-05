# Identify Hub Nodes in a Network

This function identifies hub nodes in an `igraph` network based on
degree and betweenness centrality using either z-score or quantile
thresholds. Optionally, it visualizes the classification using a scatter
plot with marginal histograms.

## Usage

``` r
find_hubs(
  graph,
  method = c("zscore", "quantile"),
  degree_threshold = 3,
  betweenness_threshold = 1,
  degree_quantile = 0.95,
  betweenness_quantile = 0.95,
  log_transform = TRUE,
  plot = TRUE,
  focus_color = "darkgreen",
  label.size = 12,
  hub_names = TRUE,
  hub_cex = 3,
  gg_extra = list(),
  weights = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph:

  An `igraph` object or a data frame containing a symbolic edge list in
  the first two columns. Additional columns are considered as edge
  attributes.

- method:

  Character. Method to identify hubs: `"zscore"` (standardized metrics)
  or `"quantile"` (empirical percentiles).

- degree_threshold:

  Numeric. Threshold for standardized degree (only used if
  `method = "zscore"`).

- betweenness_threshold:

  Numeric. Threshold for standardized betweenness (only used if
  `method = "zscore"`).

- degree_quantile:

  Numeric between 0 and 1. Quantile threshold for degree (used if
  `method = "quantile"`).

- betweenness_quantile:

  Numeric between 0 and 1. Quantile threshold for betweenness (used if
  `method = "quantile"`).

- log_transform:

  Logical. If `TRUE`, applies log-transformation to degree and
  betweenness metrics.

- plot:

  Logical. If `TRUE`, generates a plot of the degree vs. betweenness
  classification.

- focus_color:

  Character. Color to display in the focus area of the plot (hubs
  region).

- label.size:

  Numeric. Base font size for plot elements. Passed to `theme_classic`.

- hub_names:

  Logical. If `TRUE`, adds node labels to identified hubs on the plot.

- hub_cex:

  Numeric. Font size scaling factor for hub labels on the plot.

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

  A scatter plot of degree vs. betweenness with hub nodes highlighted,
  or `NULL` when `plot = FALSE`. The element is always present, so the
  return shape does not depend on the arguments. Note this is a
  `ggExtraPlot` (a grid gtable), not a `ggplot`: the scatterplot is
  passed through
  [`ggMarginal()`](https://rdrr.io/pkg/ggExtra/man/ggMarginal.html) to
  add the marginal histograms, which assembles it, so it can be printed
  but not extended with `+`. To customize it, pass layers through
  `gg_extra`; those are added before the assembly step.

- `method`:

  Description of the method and thresholds used.

- `result`:

  A `tibble` with node name, degree, betweenness, transformed metrics,
  and hub status.

- `graph`:

  The original graph with a new vertex attribute `is_hub`.

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
res <- find_hubs(g, method = "quantile", plot = FALSE)
res$method
#> [1] "Hub nodes identified by method: quantile with Degree metric threshold = 1.63293809389639 and Betweenness metric threshold = 0.213719350028092 (unweighted)"
head(res$result)
#> # A tibble: 6 × 7
#>   node  degree strength betweenness degree_metric betweenness_metric is_hub
#>   <chr>  <dbl>    <dbl>       <dbl>         <dbl>              <dbl> <lgl> 
#> 1 1         33       33      0.887           3.53             0.635  TRUE  
#> 2 2          8        8      0.237           2.20             0.213  FALSE 
#> 3 3          4        4      0.418           1.61             0.349  FALSE 
#> 4 4          7        7      0.258           2.08             0.229  TRUE  
#> 5 5          2        2      0.0253          1.10             0.0250 FALSE 
#> 6 6          2        2      0.0500          1.10             0.0488 FALSE 

# The returned graph carries an `is_hub` vertex attribute, so results chain.
table(igraph::V(res$graph)$is_hub)
#> 
#> FALSE  TRUE 
#>    77     3 
```
