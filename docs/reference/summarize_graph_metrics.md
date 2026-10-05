# Summarize Topological Properties of a Graph

Computes a comprehensive set of global topological metrics for an input
graph, including basic structure, connectivity, spectral properties, and
complexity. Supports both `igraph` objects and data frames representing
edge lists.

## Usage

``` r
summarize_graph_metrics(
  graph,
  weights = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph:

  An `igraph` object or a data frame with columns `from` and `to`
  representing an edge list.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

## Value

A one-row `data.frame`, each column a graph-level metric.

Metrics that are mathematically undefined for the input are `NaN` rather
than substituted values, as returned by the underlying igraph and ineq
functions. On degenerate graphs this is expected: an edgeless graph has
no paths (`Average_path_length`), no connected triples
(`Clustering_coefficient`), no degree variance (`Degree_assortativity`)
and a zero mean degree (`Gini_degree`).

## Details

Metrics computed:

- Number of nodes and edges

- Directed TRUE/FALSE

- Graph density

- Diameter and average path length of the largest connected component

- Clustering coefficient (transitivity)

- Degree assortativity

- Average degree and betweenness centrality

- Number of connected components and size of the largest connected
  component

- Number of single nodes

- Algebraic connectivity (second-smallest Laplacian eigenvalue)

- Degree entropy (Shannon entropy of the degree distribution)

- Gini coefficient of node degrees

- Modularity of the community structure (via Louvain algorithm)

When `weights` are supplied, the metrics switch to their weighted
definitions: `Avg_degree` becomes mean vertex strength; `Diameter`,
`Average_path_length` and `Avg_betweenness` are computed on edge
*costs*; `Clustering_coefficient` uses Barrat's weighted transitivity;
`Degree_assortativity` correlates strengths rather than degrees;
`Degree_entropy` and `Gini_degree` describe the strength distribution;
and the Laplacian behind `Algebraic_connectivity` and the Louvain run
behind `Modularity` both use edge strengths. Two columns are added,
`Is_weighted` and `Avg_strength`, so a weighted and an unweighted
summary can be row-bound.

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

- Newman, M. E. J. (2010). *Networks: An Introduction*. Oxford
  University Press.

- Estrada, E. (2012). *The Structure of Complex Networks: Theory and
  Applications*. Oxford University Press.

- Latora, V., Nicosia, V., & Russo, G. (2017). *Complex Networks:
  Principles, Methods and Applications*. Cambridge University Press.

- Louvain modularity method: Blondel, V. D., Guillaume, J. L.,
  Lambiotte, R., & Lefebvre, E. (2008). *Fast unfolding of communities
  in large networks*. J. Stat. Mech., 2008(10), P10008.

- Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A.
  (2004). *The architecture of complex weighted networks*. PNAS,
  101(11), 3747-3752.

## Examples

``` r
g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
summarize_graph_metrics(g)
#>   Nodes Edges Is_directed Is_weighted   Density Diameter Average_path_length
#> 1    60   145       FALSE       FALSE 0.0819209        6            2.768927
#>   Clustering_coefficient Degree_assortativity Avg_degree Avg_strength
#> 1             0.06280788           0.09418808   4.833333     4.833333
#>   Avg_betweenness Components Single_nodes LCC_size LCC_percent
#> 1        52.18333          1            0       60           1
#>   Algebraic_connectivity Degree_entropy Gini_degree Modularity
#> 1              0.4074096       2.991815   0.2351724  0.3614507

# Weighted: declare what the attribute means, because igraph alone would read
# it as a cost for path metrics and as a strength for modularity.
igraph::E(g)$confidence <- runif(igraph::ecount(g), 0.1, 1)
summarize_graph_metrics(g, weights = "confidence", weight_type = "strength")
#>   Nodes Edges Is_directed Is_weighted   Density Diameter Average_path_length
#> 1    60   145       FALSE        TRUE 0.0819209 19.11682            5.570991
#>   Clustering_coefficient Degree_assortativity Avg_degree Avg_strength
#> 1             0.05964938            0.1441752   4.833333     2.562503
#>   Avg_betweenness Components Single_nodes LCC_size LCC_percent
#> 1        62.91667          1            0       60           1
#>   Algebraic_connectivity Degree_entropy Gini_degree Modularity
#> 1             0.08446972       5.906891   0.2743675  0.4163293
```
