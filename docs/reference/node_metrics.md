# Compute a Table of Node-Level Centrality Metrics

Computes several per-node centrality and position metrics in one call
and returns them as a tibble, with every metric also attached to the
graph as a vertex attribute.

## Usage

``` r
node_metrics(
  graph,
  metrics = c("degree", "strength", "betweenness", "harmonic", "eigenvector", "pagerank",
    "coreness", "clustering", "constraint", "eccentricity"),
  weights = NULL,
  weight_type = c("strength", "distance"),
  normalized = TRUE,
  mode = c("all", "in", "out"),
  max_nodes = 5000,
  plot = TRUE,
  plot_type = c("correlation", "ranking"),
  top_n = 15,
  label.size = 12
)
```

## Arguments

- graph:

  An `igraph` object representing the network to analyze, or a data
  frame containing a symbolic edge list in the first two columns.
  Additional columns are considered as edge attributes.

- metrics:

  Character vector of metrics to compute. Any subset of `"degree"`,
  `"strength"`, `"betweenness"`, `"closeness"`, `"harmonic"`,
  `"eigenvector"`, `"pagerank"`, `"coreness"`, `"clustering"`,
  `"constraint"` and `"eccentricity"`. Defaults to all but `"closeness"`
  – see Details.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- normalized:

  Logical. If `TRUE` (default), metrics that igraph can normalize to
  `[0, 1]` are normalized, which makes them comparable across graphs of
  different size.

- mode:

  For directed graphs, whether to count `"all"` (default), `"in"` or
  `"out"` edges. Ignored for undirected graphs.

- max_nodes:

  Integer. Above this vertex count, the metrics whose cost is O(VE) –
  `"betweenness"`, `"closeness"`, `"harmonic"` and `"eccentricity"` –
  are skipped with a warning and returned as `NA`. Set to `Inf` to
  compute them regardless. Default is `5000`, matching the threshold
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md)
  already uses for betweenness.

- plot:

  Logical. Whether to build a diagnostic plot. Default is `TRUE`.

- plot_type:

  Either `"correlation"` (default), a Spearman correlation heatmap of
  the computed metrics, or `"ranking"`, a faceted bar chart of the
  `top_n` highest-scoring nodes per metric.

- top_n:

  Integer. Number of nodes shown per facet when `plot_type = "ranking"`.
  Default is `15`.

- label.size:

  Numeric. Base font size for plot text. Default is `12`.

## Value

A list with:

- `plot`:

  A `ggplot2` object, or `NULL` when `plot = FALSE`. The element is
  always present, so the return shape does not depend on the arguments.

- `result`:

  A tibble with one row per vertex: `node` followed by one column per
  requested metric, in the order given by `metrics`.

- `graph`:

  The input graph with each metric attached as a vertex attribute of the
  same name.

- `method`:

  A human-readable description of what was computed.

## Details

This is the node-level counterpart to
[`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md),
which describes the graph as a whole. Having both means a ranking of
nodes no longer has to be stitched together from individual igraph
calls, and because the returned graph carries the metrics as vertex
attributes, the result feeds directly into
[`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md)
(`color = "pagerank"`) and into
[`robustness_analysis()`](https://agallinat.github.io/netkit/reference/robustness_analysis.md)
(`removal_strategy = "pagerank"`), which already accepts any numeric
vertex attribute as a removal priority.

The metrics are:

- `degree`:

  Number of incident edges. Always the unweighted count.

- `strength`:

  Sum of incident edge weights; equal to `degree` when unweighted.

- `betweenness`:

  Fraction of shortest paths passing through the node. Uses edge *costs*
  when weighted.

- `closeness`:

  Reciprocal of the mean distance to all *reachable* nodes. Not
  comparable across components – see below.

- `harmonic`:

  Sum of the reciprocal distances. Unlike `closeness` this is comparable
  across components, because an unreachable pair contributes \\1/\infty
  = 0\\ rather than being excluded from the average.

- `eigenvector`:

  Leading eigenvector of the adjacency matrix: a node is central if its
  neighbors are.

- `pagerank`:

  Stationary distribution of a random surfer. Sums to 1.

- `coreness`:

  Largest *k* for which the node belongs to the k-core.

- `clustering`:

  Local transitivity: how interconnected the node's neighborhood is.
  `NaN` for nodes of degree below 2, which have no neighbor pairs.

- `constraint`:

  Burt's constraint: how much the node's connections are concentrated
  within a single group. Low constraint marks a broker.

- `eccentricity`:

  Distance to the furthest reachable node.

`"closeness"` is omitted from the defaults because it is not comparable
across components, and most real networks are disconnected. igraph
averages the distance over reachable vertices only, so a node in a
*small* isolated component scores *higher* than an equally well-placed
node in the giant component – everything in its component is nearby.
Ranking by closeness on a fragmented graph therefore puts the periphery
on top. `"harmonic"` answers the same question without that failure
mode, because unreachable pairs contribute zero rather than being
dropped from the denominator. Request `"closeness"` explicitly if the
graph is connected or you want the per-component reading.

The default plot is a correlation heatmap rather than a ranking because
the most common mistake with a table like this is to treat the metrics
as independent evidence. On many networks betweenness, eigenvector
centrality and PageRank correlate with degree above 0.9, so a node that
looks important by four measures may be important by one. The heatmap
makes that visible before the ranking is interpreted.

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

Freeman, L. C. (1978). Centrality in social networks: conceptual
clarification. *Social Networks*, 1(3), 215-239.
[doi:10.1016/0378-8733(78)90021-7](https://doi.org/10.1016/0378-8733%2878%2990021-7)

Burt, R. S. (1992). *Structural Holes: The Social Structure of
Competition*. Harvard University Press.

Marchiori, M., & Latora, V. (2000). Harmony in the small-world. *Physica
A*, 285(3-4), 539-546.
[doi:10.1016/S0378-4371(00)00311-3](https://doi.org/10.1016/S0378-4371%2800%2900311-3)

## See also

[`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md)
for the graph-level counterpart,
[`find_hubs()`](https://agallinat.github.io/netkit/reference/find_hubs.md)
for thresholded classification.

## Examples

``` r
g <- igraph::sample_pa(60, power = 1.5, directed = FALSE)
igraph::V(g)$name <- paste0("n", seq_len(igraph::vcount(g)))

res <- node_metrics(g, plot = FALSE)
head(res$result)
#> # A tibble: 6 × 11
#>   node  degree strength betweenness harmonic eigenvector pagerank coreness
#>   <chr>  <dbl>    <dbl>       <dbl>    <dbl>       <dbl>    <dbl>    <dbl>
#> 1 n1    0.186        11      0.510     0.450      0.239   0.0839         1
#> 2 n2    0.0169        1      0         0.296      0.0575  0.00898        1
#> 3 n3    0.102         6      0.571     0.451      0.391   0.0453         1
#> 4 n4    0.271        16      0.729     0.522      1       0.122          1
#> 5 n5    0.0339        2      0.0339    0.338      0.255   0.0174         1
#> 6 n6    0.0169        1      0         0.327      0.240   0.00896        1
#> # ℹ 3 more variables: clustering <dbl>, constraint <dbl>, eccentricity <dbl>

# Every metric is on the returned graph, so it chains straight onward.
igraph::vertex_attr_names(res$graph)
#>  [1] "name"         "degree"       "strength"     "betweenness"  "harmonic"    
#>  [6] "eigenvector"  "pagerank"     "coreness"     "clustering"   "constraint"  
#> [11] "eccentricity"
rob <- robustness_analysis(res$graph, removal_strategy = "pagerank",
                           steps = 10, plot = FALSE)
rob$auc$lcc_size
#> [1] 0.08637549

# The default plot shows how far the metrics actually disagree.
p <- node_metrics(g, metrics = c("degree", "betweenness", "pagerank"))$plot
```
