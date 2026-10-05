# Score an Observed Graph Against a Null Ensemble

Compares each of a graph's global metrics against the distribution that
metric takes over a null ensemble, and reports the z-score and empirical
p-value.

## Usage

``` r
metric_significance(
  graph,
  null = NULL,
  metrics = NULL,
  n = 100,
  model = c("rewire", "configuration", "erdos_renyi"),
  weights = NULL,
  weight_type = c("strength", "distance"),
  seed = NULL,
  plot = TRUE,
  label.size = 12
)
```

## Arguments

- graph:

  An `igraph` object or a data frame edge list, as elsewhere in netkit.

- null:

  A list of null graphs, as produced by
  [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md).
  If `NULL` (default), an ensemble is generated internally with `n`
  graphs under `model`.

- metrics:

  Character vector of metric names to test, matching columns of
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md).
  If `NULL` (default), every numeric metric that varies across the
  ensemble is used.

- n:

  Integer. Number of null graphs to generate when `null` is `NULL`.
  Default is `100`. Named `n` rather than `n_null` so that it matches
  [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
  and so that `n =` is an exact argument match – with a formal called
  `n_null` alongside `null`, `n =` is ambiguous and errors.

- model:

  Passed to
  [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
  when `null` is `NULL`.

- weights, weight_type:

  Passed to
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md).
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- seed:

  Integer or `NULL`. Seed for reproducibility.

- plot:

  Logical. Whether to build a diagnostic plot. Default is `TRUE`.

- label.size:

  Numeric. Base font size for plot text. Default is `12`.

## Value

A list with:

- `plot`:

  A faceted `ggplot2` object showing each metric's null distribution
  with the observed value marked, or `NULL` when `plot = FALSE`. The
  element is always present.

- `result`:

  A tibble with one row per metric: `metric`, `observed`, `null_mean`,
  `null_sd`, `z`, `p_empirical`, `ci_lower` and `ci_upper` (the 2.5th
  and 97.5th percentiles of the null).

- `graph`:

  The input graph, unchanged, so the call still chains.

- `method`:

  A human-readable description of the test performed.

## Details

`p_empirical` is two-sided and uses the `(r + 1) / (n + 1)` convention,
where `r` counts null values at least as extreme as the observed one. It
is therefore never exactly zero: with `n` nulls the smallest attainable
p-value is `1 / (n + 1)`, so testing against 20 nulls cannot produce
evidence at `p < 0.05` however large the effect. Choose `n` accordingly.

`z` is `NA` when the null distribution has zero variance – which happens
legitimately for metrics the null model holds fixed, such as `Nodes`,
`Edges` and `Density` under a degree-preserving model. Those metrics are
excluded from the automatic `metrics` selection for that reason, but are
reported as `NA` rather than dropped if you request them explicitly.

## See also

[`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
for the ensembles,
[`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md)
for the metrics themselves.

## Examples

``` r
g <- igraph::sample_pa(60, power = 1.5, directed = FALSE)

# `n` is small here to keep the example fast. Note the floor this puts
# on the attainable p-value: 1 / (10 + 1).
res <- metric_significance(g, metrics = c("Clustering_coefficient",
                                          "Modularity"),
                           n = 10, seed = 1, plot = FALSE)
res$result
#> # A tibble: 2 × 8
#>   metric          observed null_mean null_sd     z p_empirical ci_lower ci_upper
#>   <chr>              <dbl>     <dbl>   <dbl> <dbl>       <dbl>    <dbl>    <dbl>
#> 1 Clustering_coe…    0         0.180  0.0920 -1.96      0.0909   0.0300    0.321
#> 2 Modularity         0.663     0.625  0.0121  3.13      0.0909   0.604     0.637
```
