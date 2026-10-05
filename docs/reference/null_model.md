# Generate Null-Model Graphs for Significance Testing

Produces an ensemble of random graphs matched to an observed graph in
some respect, so that an observed metric can be compared against what
that respect alone would predict.

## Usage

``` r
null_model(
  graph,
  model = c("rewire", "configuration", "erdos_renyi"),
  n = 100,
  shuffle_weights = TRUE,
  seed = NULL,
  niter = NULL
)
```

## Arguments

- graph:

  An `igraph` object representing the network to analyze, or a data
  frame containing a symbolic edge list in the first two columns.
  Additional columns are considered as edge attributes.

- model:

  Which aspect of the observed graph to preserve:

  `"rewire"`

  :   Degree-preserving edge rewiring (the default). Preserves the
      degree sequence exactly while staying close to the observed graph.

  `"configuration"`

  :   Draws a fresh graph with the same degree sequence via
      [`igraph::sample_degseq()`](https://r.igraph.org/reference/sample_degseq.html).
      Preserves the degree sequence exactly but is otherwise
      unconstrained.

  `"erdos_renyi"`

  :   Matches only the vertex and edge counts. The weakest null,
      included as the deliberate contrast: comparing against it tells
      you whether an effect needs anything beyond density.

- n:

  Integer. Number of null graphs to generate. Default is `100`.

- shuffle_weights:

  Logical. If `TRUE` (default) and the observed graph carries a `weight`
  edge attribute, the observed weights are permuted across the null
  graph's edges, so the weight *distribution* is preserved while its
  placement is randomized. If `FALSE`, the null graphs carry no weights
  at all – including under `model = "rewire"`, which would otherwise
  inherit the observed attribute bound to edges that have since been
  rewired.

- seed:

  Integer or `NULL`. Seed for reproducibility. When `NULL` (default) the
  caller's random stream is used and left alone.

- niter:

  Integer or `NULL`. Number of rewiring trials for `model = "rewire"`.
  When `NULL` (default), `10 * ecount(graph)`, which is the usual rule
  of thumb for mixing.

## Value

A list of `n` `igraph` objects, with class `c("netkit_null", "list")`
and the model recorded in its `"model"` attribute.

## Details

This is what turns
[`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md)
from descriptive into inferential. On its own, a modularity of 0.42 says
nothing: random graphs with the same degree sequence routinely reach
0.3. Pass an ensemble from this function to
[`metric_significance()`](https://agallinat.github.io/netkit/reference/metric_significance.md)
to get the z-score and empirical p-value.

`"rewire"` and `"configuration"` both preserve the degree sequence,
which is almost always the right null for a topological claim: nearly
every network metric is partly determined by the degree sequence, so a
null that does not hold it fixed mostly rediscovers that the graph is
heavy-tailed.

They differ in how far they travel from the observed graph. `"rewire"`
makes local double-edge swaps, so the ensemble is centered on the
observed graph and is the more conservative choice. `"configuration"`
resamples from scratch.

`"configuration"` is attempted with `method = "vl"`, which produces
simple connected graphs, and falls back to `"configuration.simple"` with
a warning where that is not possible – `"vl"` requires a connected
realization of the degree sequence to exist, which fails for example
when the graph has isolated vertices.

## References

Maslov, S., & Sneppen, K. (2002). Specificity and stability in topology
of protein networks. *Science*, 296(5569), 910-913.
[doi:10.1126/science.1065103](https://doi.org/10.1126/science.1065103)

Newman, M. E. J., Strogatz, S. H., & Watts, D. J. (2001). Random graphs
with arbitrary degree distributions and their applications. *Physical
Review E*, 64(2), 026118.
[doi:10.1103/PhysRevE.64.026118](https://doi.org/10.1103/PhysRevE.64.026118)

## See also

[`metric_significance()`](https://agallinat.github.io/netkit/reference/metric_significance.md)
to score an observed graph against the ensemble,
[`small_worldness()`](https://agallinat.github.io/netkit/reference/small_worldness.md)
for the classic application.

## Examples

``` r
g <- igraph::sample_pa(60, power = 1.5, directed = FALSE)

# `n` is small here to keep the example fast; use 100 or more in practice.
nulls <- null_model(g, model = "rewire", n = 10, seed = 1)
length(nulls)
#> [1] 10

# The degree sequence is preserved exactly, which is the point.
identical(sort(igraph::degree(nulls[[1]])), sort(igraph::degree(g)))
#> [1] TRUE
```
