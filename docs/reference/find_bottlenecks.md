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
  gg_extra = list()
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

## Examples

``` r
g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
res <- find_bottlenecks(g, method = "quantile", plot = FALSE)
res$method
#> [1] "Bottlenecks identified by method: quantile with Degree metric threshold = 0.693147180559945 and Betweenness metric threshold = 0.00625032555135432"
head(res$result)
#> # A tibble: 6 × 6
#>   node  degree betweenness degree_metric betweenness_metric is_bottleneck
#>   <chr>  <dbl>       <dbl>         <dbl>              <dbl> <lgl>        
#> 1 1          4      0.235          1.61              0.211  FALSE        
#> 2 2         46      0.967          3.85              0.676  FALSE        
#> 3 3          2      0.0500         1.10              0.0488 FALSE        
#> 4 4          1      0              0.693             0      FALSE        
#> 5 5          3      0.0503         1.39              0.0491 FALSE        
#> 6 6          2      0.0500         1.10              0.0488 FALSE        
```
