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
  label.size = 12
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

## Value

A list with:

- `plot`:

  A `ggplot2` object overlaying both CCDF curves.

- `global_topology`:

  A data frame with one row per input graph, as produced by
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md).

- `similarity`:

  A one-row data frame with `jaccard_similarity` (Jaccard index of the
  edge sets), `node_overlap` (fraction of shared nodes) and
  `edge_overlap` (fraction of shared edges).

- `ks_test`:

  The Kolmogorov-Smirnov test comparing the two degree distributions, as
  returned by
  [`stats::ks.test()`](https://rdrr.io/r/stats/ks.test.html).

## Examples

``` r
g1 <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
g2 <- igraph::sample_gnp(80, 0.05, directed = FALSE)
res <- compare_networks(g1, g2)
res$global_topology
#>   Nodes Edges Is_directed    Density Diameter Average_path_length
#> 1    80    79       FALSE 0.02500000        7             3.21519
#> 2    80   179       FALSE 0.05664557        7             3.09088
#>   Clustering_coefficient Degree_assortativity Avg_degree Avg_betweenness
#> 1             0.00000000          -0.38504587      1.975          87.500
#> 2             0.07939683           0.02884186      4.475          80.525
#>   Components Single_nodes LCC_size LCC_percent Algebraic_connectivity
#> 1          1            0       80      1.0000           7.872572e-02
#> 2          2            1       79      0.9875           6.990798e-18
#>   Degree_entropy Gini_degree Modularity
#> 1       1.399238   0.4359177  0.6154462
#> 2       3.052382   0.2655028  0.4223495
res$ks_test
#> 
#>  Exact two-sample Kolmogorov-Smirnov test
#> 
#> data:  deg1 and deg2
#> D = 0.7, p-value < 2.2e-16
#> alternative hypothesis: two-sided
#> 
```
