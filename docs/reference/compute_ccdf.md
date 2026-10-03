# Compute the Complementary Cumulative Distribution Function (CCDF) of Node Degrees

Computes the CCDF of node degrees for a given igraph object. The CCDF is
useful for visualizing degree distributions, particularly on log-log
plots, to identify power-law or heavy-tailed behaviors. Internal helper
function for
[`plot_CCDF()`](https://agallinat.github.io/netkit/reference/plot_CCDF.md)
and
[`compare_networks()`](https://agallinat.github.io/netkit/reference/compare_networks.md).

## Usage

``` r
compute_ccdf(graph, mode = c("all", "in", "out"), remove_singles = FALSE)
```

## Arguments

- graph:

  An igraph object representing the graph.

- mode:

  Character string indicating which degree type to compute. Options are
  `"all"` (default), `"in"`, or `"out"`. Only relevant for directed
  graphs.

- remove_singles:

  Logical. If `TRUE`, nodes with degree zero will be removed before
  computing the degree distribution.

## Value

A data frame with two columns:

- degree:

  Integer node degree values.

- ccdf:

  Complementary cumulative distribution values (P(X ≥ x)).
