# Extract an Interpretable Subnetwork Around a Set of Nodes

Given a set of nodes of interest – a gene list, a set of disease genes,
a group of drug targets – returns a connected, readable subnetwork
around them. Five strategies are offered, from the plain induced
subgraph through to an approximate Steiner tree and a diffusion-based
expansion.

## Usage

``` r
extract_subnetwork(
  graph,
  nodes,
  method = c("induced", "neighbors", "shortest_paths", "steiner", "diffusion"),
  order = 1,
  max_degree = NULL,
  top_n = 100,
  diffusion_method = c("rwr", "laplacian", "heat"),
  largest_component = FALSE,
  weights = NULL,
  weight_type = c("strength", "distance"),
  plot = TRUE,
  ...
)
```

## Arguments

- graph:

  An `igraph` object representing the network to analyze, or a data
  frame containing a symbolic edge list in the first two columns.
  Additional columns are considered as edge attributes.

- nodes:

  Character vector of vertex names to build the subnetwork around (the
  "seeds", or in the Steiner case the "terminals").

- method:

  How to choose the nodes to keep:

  `"induced"`

  :   Only `nodes` themselves, with whatever edges run between them. The
      baseline every other method should be compared against.

  `"neighbors"`

  :   `nodes` plus their `order`-step neighborhood.

  `"shortest_paths"`

  :   The union of all shortest paths between every pair of `nodes`. The
      classic "connect my gene list".

  `"steiner"`

  :   An approximate minimum Steiner tree over `nodes`: the smallest
      tree that connects them all. Usually the most readable figure,
      because it is a tree rather than a union of overlapping paths.

  `"diffusion"`

  :   Propagate from `nodes` with
      [`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md)
      and keep the `top_n` highest-scoring vertices.

- order:

  Integer. Neighborhood radius when `method = "neighbors"`. Default is
  `1`.

- max_degree:

  Integer or `NULL`. When `method = "neighbors"`, exclude neighbors
  whose degree exceeds this. Without it, first-neighbor expansion on a
  hub-dominated network returns most of the network; see Details.

- top_n:

  Integer. Number of vertices to keep when `method = "diffusion"`.
  Default is `100`.

- diffusion_method:

  Passed to
  [`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md)
  as `method` when `method = "diffusion"`. Default is `"rwr"`.

- largest_component:

  Logical. If `TRUE`, return only the largest connected component of the
  result. Default is `FALSE`, which keeps a forest when the seeds span
  several components.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  Path lengths and the Steiner tree use edge *costs*; diffusion uses
  *strengths*. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- plot:

  Logical. If `TRUE` (default), draws the extracted subnetwork with
  [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md),
  seeds highlighted.

- ...:

  Additional arguments passed to
  [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md).

## Value

A list with:

- `result`:

  A tibble with one row per retained vertex: `node`, `reason` (why it
  was kept – `"seed"`, `"neighbor"`, `"on_path"`, `"steiner"` or
  `"diffused"`), `is_seed`, and `score` (the diffusion score, or `NA`
  for the other methods).

- `graph`:

  The extracted subgraph, with `is_seed` and `reason` attached as vertex
  attributes. For `method = "steiner"` this is the Steiner tree itself –
  exactly the tree's edges. For every other method it is the subgraph
  *induced* on the selected vertices, so it carries all edges that run
  between them, which is more context but may contain cycles.

- `method`:

  A human-readable description of what was done.

Note there is no `plot` element: like
[`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md),
[`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md)
and
[`highlight_nodes()`](https://agallinat.github.io/netkit/reference/highlight_nodes.md),
this function renders through base graphics via `plot.igraph` and so has
no plot object to hand back.

## Details

`"shortest_paths"` and `"steiner"` both answer "how are these nodes
connected to each other", but differently. The union of shortest paths
keeps every equally short route, so it grows quickly and can be dense;
the Steiner tree keeps one connecting structure of minimum total cost,
so it is always a tree with `vcount - 1` edges and is far easier to
read. Use the union when you care about redundancy of connection, the
tree when you want a figure.

The Steiner tree problem is NP-hard. This implementation uses the
Kou-Markowsky-Berman heuristic: build the metric closure over the
terminals (their pairwise shortest-path distances), take a minimum
spanning tree of that closure, expand each closure edge back into the
actual path it stands for, and prune non-terminal leaves. The result is
guaranteed within a factor of `2 - 2/|terminals|` of the true optimum,
which in practice is close.

Because the problem is only solved approximately, the tree returned is
one of possibly several near-minimal trees, and which one depends on how
ties between equal-cost paths are broken. That tie-breaking is not the
same in the weighted and unweighted code paths – igraph uses
breadth-first search when no weights are given and Dijkstra when they
are – so passing `weights` whose values happen to be all equal can
return a different tree of slightly different total cost than omitting
them. Both satisfy the approximation bound; neither is "the" Steiner
tree. Treat the specific vertex set as one valid answer rather than a
canonical one.

`max_degree` exists because first-neighbor expansion is dominated by
hubs. In a protein interaction network a handful of promiscuous proteins
are adjacent to a large fraction of the graph, so the one-step
neighborhood of almost any gene list is most of the network. Capping
degree removes them and leaves a subnetwork whose edges carry
information.

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

Kou, L., Markowsky, G., & Berman, L. (1981). A fast algorithm for
Steiner trees. *Acta Informatica*, 15(2), 141-145.
[doi:10.1007/BF00288961](https://doi.org/10.1007/BF00288961)

## See also

[`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md)
for the scores behind `method = "diffusion"`,
[`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md)
for unsupervised community structure.

## Examples

``` r
g <- igraph::sample_pa(80, power = 1.5, directed = FALSE)
igraph::V(g)$name <- paste0("n", seq_len(igraph::vcount(g)))
seeds <- c("n5", "n20", "n40", "n60")

# The minimal tree connecting the seeds: always vcount - 1 edges.
st <- extract_subnetwork(g, seeds, method = "steiner", plot = FALSE)
igraph::vcount(st$graph)
#> [1] 5
st$result
#> # A tibble: 5 × 4
#>   node  reason  is_seed score
#>   <chr> <chr>   <lgl>   <dbl>
#> 1 n2    steiner FALSE      NA
#> 2 n5    seed    TRUE       NA
#> 3 n20   seed    TRUE       NA
#> 4 n40   seed    TRUE       NA
#> 5 n60   seed    TRUE       NA

# Every shortest route between them, which is a superset.
sp <- extract_subnetwork(g, seeds, method = "shortest_paths", plot = FALSE)
igraph::vcount(sp$graph) >= igraph::vcount(st$graph)
#> [1] TRUE

# `top_n` is small here to keep the example fast on an 80-node graph.
df <- extract_subnetwork(g, seeds, method = "diffusion", top_n = 15,
                         plot = FALSE)
head(df$result)
#> # A tibble: 6 × 4
#>   node  reason   is_seed  score
#>   <chr> <chr>    <lgl>    <dbl>
#> 1 n2    diffused FALSE   0.0231
#> 2 n4    diffused FALSE   0.0162
#> 3 n5    seed     TRUE    0.554 
#> 4 n12   diffused FALSE   0.0162
#> 5 n15   diffused FALSE   0.0162
#> 6 n17   diffused FALSE   0.388 
```
