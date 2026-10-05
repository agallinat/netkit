# Detect and Visualize Network Modules (Communities)

Identifies modules (communities) in a network using a variety of
community detection algorithms from the igraph package. Optionally
filters out small modules, visualizes the detected modules, and returns
induced subgraphs for each module.

## Usage

``` r
find_modules(
  graph,
  method = "louvain",
  min_size = 3,
  no.of.communities = NULL,
  return_subgraphs = FALSE,
  plot = TRUE,
  label = FALSE,
  weights = NULL,
  weight_type = c("strength", "distance"),
  ...
)
```

## Arguments

- graph:

  An `igraph` object representing the network to analyze or a data frame
  containing a symbolic edge list in the first two columns. Additional
  columns are considered as edge attributes.

- method:

  Character. Community detection method. Options include: `"louvain"`,
  `"walktrap"`, `"infomap"`, `"edge_betweenness"`,
  `"fluid_communities"`, `"fast_greedy"`, `"leading_eigen"`, `"leiden"`,
  and `"spinglass"`.

- min_size:

  Integer. Minimum number of nodes required to retain a module. Modules
  smaller than this size are discarded. Default is 3.

- no.of.communities:

  Integer. Required only when `method = "fluid_communities"`. Specifies
  the number of communities to find.

- return_subgraphs:

  Logical. If `TRUE`, returns a list of induced subgraphs for each
  detected module.

- plot:

  Logical. If `TRUE`, generates a network plot colored by module.

- label:

  Logical. If `TRUE`, displays node labels in the plot.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  Community detection reads a weight as a *strength*, except
  `method = "edge_betweenness"`, which needs costs and is given them.
  `method = "fluid_communities"` cannot use weights and warns. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- ...:

  Additional parameters passed to the
  [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md)
  function for customizing the plot.

## Value

A list with the following components:

- `result`:

  A tibble mapping each node to its module assignment.

- `module_table`:

  Deprecated alias for `result`, kept for backward compatibility.

- `n_modules`:

  The number of modules that meet the `min_size` threshold.

- `subgraphs`:

  A named list of subgraphs for each module (only if
  `return_subgraphs = TRUE`).

- `method`:

  The community detection method used.

- `graph`:

  The input graph with assigned 'module' and 'color' as vertex
  attributes, if `plot = TRUE`.

If `plot = TRUE`, a network plot is displayed with nodes colored by
module.

## Details

This function is a wrapper around several igraph community detection
algorithms, including Louvain (`cluster_louvain()`), Walktrap, Infomap,
Fast Greedy, and others. It simplifies their application and offers
optional filtering, visualization via
[`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md),
and module subgraph extraction.

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

## References

Csardi G, Nepusz T. The igraph software package for complex network
research. InterJournal, Complex Systems. 2006;1695. <https://igraph.org>

## Examples

``` r
g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
res <- find_modules(g, method = "louvain", plot = FALSE)
res$n_modules
#> [1] 8
head(res$module_table)
#> # A tibble: 6 × 2
#>   node  module
#>   <chr>  <int>
#> 1 1          1
#> 2 2          2
#> 3 5          5
#> 4 6          6
#> 5 8          2
#> 6 9          6
```
