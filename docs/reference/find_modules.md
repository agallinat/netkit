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

- ...:

  Additional parameters passed to the
  [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md)
  function for customizing the plot.

## Value

A list with the following components:

- `module_table`:

  A tibble mapping each node to its module assignment.

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

## References

Csardi G, Nepusz T. The igraph software package for complex network
research. InterJournal, Complex Systems. 2006;1695. <https://igraph.org>

## Examples

``` r
g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
res <- find_modules(g, method = "louvain", plot = FALSE)
res$n_modules
#> [1] 6
head(res$module_table)
#> # A tibble: 6 × 2
#>   node  module
#>   <chr>  <int>
#> 1 1          1
#> 2 2          2
#> 3 3          3
#> 4 4          4
#> 5 5          3
#> 6 6          2
```
