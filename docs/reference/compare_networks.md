# Compare Two Networks

This function compares two networks using summary metrics, degree
distributions, and topological similarity measures. It overlays the
complementary cumulative frequency distributions (CCDFs) of degree and
returns a combined report.

## Usage

``` r
compare_networks(
  graph1,
  graph2,
  remove_singles = FALSE,
  show_PL = TRUE,
  PL_exponents = c(2, 3),
  colors = c("#e41a1c", "#000831", "#9c52f2", "#b8b8ff"),
  label.size = 12,
  weights1 = NULL,
  weights2 = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph1:

  An igraph object or data.frame (edge list).

- graph2:

  An igraph object or data.frame (edge list).

- remove_singles:

  Logical; remove single nodes before analysis.

- show_PL:

  Logical; whether to fit and display power law exponents.

- PL_exponents:

  Vector; power-law slopes to show.

- colors:

  Optional vector of colors for the CCDF plot.

- label.size:

  Labels' size in the CCDF plot.

- weights1, weights2:

  Optional edge weights for each graph: `NULL` (default) to ignore them,
  the name of an edge attribute, or a numeric vector. Supplied
  separately because the two graphs need not carry the same attribute.
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`, applied to both graphs.
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

## Value

A list with:

- `plot`:

  A `ggplot2` object overlaying both CCDF curves.

- `global_topology`:

  A data frame with one row per input graph, as produced by
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md).

- `similarity`:

  A one-row data frame with three distinct set statistics:
  `jaccard_similarity`, the Jaccard index of the edge sets (shared edges
  over their union); `node_overlap`, the Jaccard index of the vertex
  sets; and `edge_overlap`, the overlap coefficient of the edge sets
  (shared edges over the *smaller* of the two edge sets). The last is
  the informative companion to Jaccard when the two networks differ
  greatly in size – a small network nested inside a large one scores
  near 1 on overlap and near 0 on Jaccard. Each is `NaN` where its
  denominator is empty.

- `ks_test`:

  The Kolmogorov-Smirnov test comparing the two degree distributions, as
  returned by
  [`stats::ks.test()`](https://rdrr.io/r/stats/ks.test.html).

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

## Examples

``` r
g1 <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
g2 <- igraph::sample_gnp(80, 0.05, directed = FALSE)
res <- compare_networks(g1, g2)
res$global_topology
#>   Nodes Edges Is_directed Is_weighted    Density Diameter Average_path_length
#> 1    80    79       FALSE       FALSE 0.02500000        7             3.21519
#> 2    80   179       FALSE       FALSE 0.05664557        7             3.09088
#>   Clustering_coefficient Degree_assortativity Avg_degree Avg_strength
#> 1             0.00000000          -0.38504587      1.975        1.975
#> 2             0.07939683           0.02884186      4.475        4.475
#>   Avg_betweenness Components Single_nodes LCC_size LCC_percent
#> 1          87.500          1            0       80      1.0000
#> 2          80.525          2            1       79      0.9875
#>   Algebraic_connectivity Degree_entropy Gini_degree Modularity
#> 1           7.872572e-02       1.399238   0.4359177  0.6154462
#> 2           6.990798e-18       3.052382   0.2655028  0.4223495
res$ks_test
#> 
#>  Exact two-sample Kolmogorov-Smirnov test
#> 
#> data:  deg1 and deg2
#> D = 0.7, p-value < 2.2e-16
#> alternative hypothesis: two-sided
#> 
```
