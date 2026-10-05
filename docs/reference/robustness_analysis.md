# Network Robustness Analysis via Node Removal Simulation

Simulates the removal of nodes from a network using various strategies
and evaluates how the structure degrades using selected robustness
metrics. Useful for assessing the vulnerability or resilience of a
graph.

## Usage

``` r
robustness_analysis(
  graph,
  removal_strategy = c("random", "degree", "betweenness", "strength"),
  steps = 50,
  metrics = c("lcc_size", "efficiency", "n_components"),
  n_reps = 50,
  plot = TRUE,
  seed = NULL,
  weights = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph:

  An `igraph` object representing the network to plot or a data frame
  containing a symbolic edge list in the first two columns. Additional
  columns are considered as edge attributes. Must be undirected;
  directed graphs will be converted.

- removal_strategy:

  Character. Strategy used for node removal. Options are: `"random"`,
  `"degree"`, `"betweenness"`, `"strength"`, or the name of a numeric
  vertex attribute. Custom attributes are interpreted as priority scores
  (higher = removed first). `"strength"` requires `weights` and errors
  without them.

- steps:

  Integer. Number of removal steps (default: 50).

- metrics:

  Character vector. Structural metrics to compute at each step. Options
  include: `"lcc_size"`, `"efficiency"`, and `"n_components"`.

- n_reps:

  Integer. Number of simulation repetitions (only relevant if
  `removal_strategy = "random"`).

- plot:

  Logical. If `TRUE`, a robustness plot is generated.

- seed:

  Integer or NULL. Random seed for reproducibility.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  When supplied, global efficiency is computed over weighted path
  lengths, `removal_strategy = "betweenness"` uses cost-based
  betweenness, and `removal_strategy = "strength"` becomes available.
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

## Value

A list with:

- `plot`:

  A `ggplot2` object showing the evolution of the selected metrics as
  nodes are progressively removed, or `NULL` when `plot = FALSE`. The
  element is always present, so the return shape does not depend on the
  arguments.

- `result`:

  A summarized data frame (mean and SD) if `n_reps > 1`, otherwise raw
  results.

- `all_results`:

  A data frame with simulation results across all steps and repetitions.

- `summary`:

  Deprecated alias for `result`, kept for backward compatibility.

- `auc`:

  Named list of AUC (area under the curve) values for each selected
  metric.

## Details

This function builds on classic approaches in network science for
evaluating structural robustness, simulating progressive node removal
and quantifying the degradation of key topological features.

For deterministic strategies (`"degree"`, `"betweenness"`, or custom
attributes), nodes are removed in a fixed priority order. For the
`"random"` strategy, the process is repeated `n_reps` times, and the
results are aggregated.

The available metrics are:

- **Largest Connected Component**: size of the largest remaining
  component (`lcc_size`).

- **Global Efficiency**: average inverse shortest path length among all
  pairs (`efficiency`).

- **Number of Components**: total number of disconnected components
  (`n_components`).

Additionally, Area Under the Curve (AUC) is calculated for each metric,
providing a scalar summary of robustness. A higher AUC indicates greater
resilience (i.e., slower degradation).

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

Albert R, Jeong H, Barabási AL. Error and attack tolerance of complex
networks. Nature. 2000;406(6794):378–382.
[doi:10.1038/35019019](https://doi.org/10.1038/35019019)

## Examples

``` r
g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)

# `steps` is reduced from its default of 50 to keep the example fast.
res <- robustness_analysis(g, removal_strategy = "degree", steps = 10,
                           metrics = c("lcc_size", "n_components"),
                           plot = FALSE, seed = 1)
res$auc
#> $lcc_size
#> [1] 0.1132208
#> 
#> $n_components
#> [1] 0.5945002
#> 
head(res$summary)
#> # A tibble: 6 × 5
#>     rep removed removed_frac lcc_size n_components
#>   <int>   <int>        <dbl>    <dbl>        <dbl>
#> 1     1       1       0.0127       18           25
#> 2     1       8       0.101         4           55
#> 3     1      16       0.203         2           58
#> 4     1      24       0.304         1           56
#> 5     1      32       0.405         1           48
#> 6     1      40       0.506         1           40
```
