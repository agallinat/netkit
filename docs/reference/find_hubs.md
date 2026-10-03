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
  gg_extra = list()
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

## Examples

``` r
g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
res <- find_hubs(g, method = "quantile", plot = FALSE)
res$method
#> [1] "Hub nodes identified by method: quantile with Degree metric threshold = 1.79175946922805 and Betweenness metric threshold = 0.117872475763818"
head(res$result)
#> # A tibble: 6 × 6
#>   node  degree betweenness degree_metric betweenness_metric is_hub
#>   <chr>  <dbl>       <dbl>         <dbl>              <dbl> <lgl> 
#> 1 1         28      0.951          3.37              0.668  TRUE  
#> 2 2          5      0.123          1.79              0.116  FALSE 
#> 3 3          3      0.0987         1.39              0.0941 FALSE 
#> 4 4          4      0.0750         1.61              0.0723 FALSE 
#> 5 5          2      0.0253         1.10              0.0250 FALSE 
#> 6 6          1      0              0.693             0      FALSE 

# The returned graph carries an `is_hub` vertex attribute, so results chain.
table(igraph::V(res$graph)$is_hub)
#> 
#> FALSE  TRUE 
#>    78     2 
```
