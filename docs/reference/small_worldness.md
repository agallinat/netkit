# Small-World Coefficients

Computes the small-world coefficient sigma, which compares a graph's
clustering and path length against a degree-matched random ensemble.

## Usage

``` r
small_worldness(
  graph,
  n_null = 100,
  model = c("rewire", "configuration", "erdos_renyi"),
  weights = NULL,
  weight_type = c("strength", "distance"),
  seed = NULL
)
```

## Arguments

- graph:

  An `igraph` object or a data frame edge list.

- n_null:

  Integer. Number of null graphs. Default is `100`.

- model:

  Passed to
  [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md).
  Default is `"rewire"`.

- weights, weight_type:

  Passed to
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md).
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- seed:

  Integer or `NULL`. Seed for reproducibility.

## Value

A list with:

- `result`:

  A one-row tibble: `sigma`, `C`, `C_rand`, `L`, `L_rand` and `n_null`.

- `graph`:

  The input graph, unchanged.

- `method`:

  A human-readable description.

## Details

A network is "small-world" when it is much more clustered than a random
graph with the same degree sequence while having a comparable average
path length. `sigma` expresses that as a single ratio:
`(C/C_rand) / (L/L_rand)`, where values appreciably above 1 indicate
small-world organization.

Only `sigma` is reported. The companion coefficient `omega` of Telesford
et al. (2011) additionally requires a *lattice* reference, and igraph
provides no degree-preserving latticization; implementing one
approximately would make `omega` quietly dependent on how well that
approximation worked, so it is omitted rather than shipped unreliable.

`sigma` is known to grow with network size, so it is not comparable
across graphs of different size – use it to ask whether one graph is
small-world, not which of two is more so.

## References

Humphries, M. D., & Gurney, K. (2008). Network "small-world-ness": a
quantitative method for determining canonical network equivalence. *PLoS
ONE*, 3(4), e0002051.
[doi:10.1371/journal.pone.0002051](https://doi.org/10.1371/journal.pone.0002051)

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
"small-world" networks. *Nature*, 393(6684), 440-442.
[doi:10.1038/30918](https://doi.org/10.1038/30918)

## See also

[`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md),
[`metric_significance()`](https://agallinat.github.io/netkit/reference/metric_significance.md)

## Examples

``` r
g <- igraph::sample_smallworld(1, 60, 4, 0.05)

# `n_null` is small here to keep the example fast.
sw <- small_worldness(g, n_null = 10, seed = 1)
sw$result
#> # A tibble: 1 × 6
#>   sigma     C C_rand     L L_rand n_null
#>   <dbl> <dbl>  <dbl> <dbl>  <dbl>  <int>
#> 1  3.40 0.470  0.114  2.61   2.15     10
```
