# Network Robustness Analysis via Node Removal Simulation

Simulates the removal of nodes from a network using various strategies
and evaluates how the structure degrades using selected robustness
metrics. Useful for assessing the vulnerability or resilience of a
graph.

## Usage

``` r
robustness_analysis(
  graph,
  removal_strategy = c("random", "degree", "betweenness"),
  steps = 50,
  metrics = c("lcc_size", "efficiency", "n_components"),
  n_reps = 50,
  plot = TRUE,
  seed = NULL
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
  `"degree"`, `"betweenness"`, or the name of a numeric vertex
  attribute. Custom attributes are interpreted as priority scores
  (higher = removed first).

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

## Value

A list with:

- `plot`:

  A `ggplot2` object showing the evolution of the selected metrics as
  nodes are progressively removed, or `NULL` when `plot = FALSE`. The
  element is always present, so the return shape does not depend on the
  arguments.

- `all_results`:

  A data frame with simulation results across all steps and repetitions.

- `summary`:

  A summarized data frame (mean and SD) if `n_reps > 1`, otherwise raw
  results.

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
#> [1] 0.1099267
#> 
#> $n_components
#> [1] 0.5946149
#> 
head(res$summary)
#> # A tibble: 6 × 5
#>     rep removed removed_frac lcc_size n_components
#>   <int>   <int>        <dbl>    <dbl>        <dbl>
#> 1     1       1       0.0127       19           30
#> 2     1       8       0.101         3           59
#> 3     1      16       0.203         3           58
#> 4     1      24       0.304         1           56
#> 5     1      32       0.405         1           48
#> 6     1      40       0.506         1           40
```
