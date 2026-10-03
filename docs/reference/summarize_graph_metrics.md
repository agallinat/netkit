# Summarize Topological Properties of a Graph

Computes a comprehensive set of global topological metrics for an input
graph, including basic structure, connectivity, spectral properties, and
complexity. Supports both `igraph` objects and data frames representing
edge lists.

## Usage

``` r
summarize_graph_metrics(graph)
```

## Arguments

- graph:

  An `igraph` object or a data frame with columns `from` and `to`
  representing an edge list.

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

## Examples

``` r
g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
summarize_graph_metrics(g)
#>   Nodes Edges Is_directed   Density Diameter Average_path_length
#> 1    60   131       FALSE 0.0740113        6            2.983616
#>   Clustering_coefficient Degree_assortativity Avg_degree Avg_betweenness
#> 1             0.07453476           0.06304963   4.366667        58.51667
#>   Components Single_nodes LCC_size LCC_percent Algebraic_connectivity
#> 1          1            0       60           1              0.4683649
#>   Degree_entropy Gini_degree Modularity
#> 1       2.756289   0.2279898  0.4186528
```
