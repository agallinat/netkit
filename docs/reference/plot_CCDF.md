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
  label.size = 12
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

## Value

A `ggplot2` object showing the CCDF of node degrees on a log-log scale.

## Examples

``` r
g <- igraph::sample_pa(200, power = 1.5, directed = FALSE)
plot_CCDF(g)


# Compare against reference power-law slopes.
plot_CCDF(g, PL_exponents = c(2, 2.5, 3))

```
