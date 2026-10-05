# Greedy Seed Node Selection to Maximize Diffusion Toward Target Nodes

This function implements a greedy algorithm to select a set of `k` seed
nodes from a candidate list such that the resulting diffusion signal on
a specified set of target nodes is maximized. It supports multiple
diffusion models (Laplacian, Heat Kernel, Random Walk with Restart).

## Usage

``` r
greedy_seed_selection(
  graph,
  target_nodes,
  candidate_nodes = NULL,
  method = c("laplacian", "heat", "rwr"),
  alpha = 0.7,
  t = 1,
  restart_prob = 0.3,
  k = 5,
  normalize = TRUE,
  plot = TRUE,
  weights = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph:

  An `igraph` object representing the network to analyze or a data frame
  containing a symbolic edge list in the first two columns. Additional
  columns are considered as edge attributes.

- target_nodes:

  A character vector of node names to prioritize for receiving the
  diffusion signal.

- candidate_nodes:

  Optional character vector of eligible nodes to consider as seeds. If
  NULL (default), all non-target nodes are used.

- method:

  Diffusion method to use. One of `"laplacian"`, `"heat"`, or `"rwr"`.

- alpha:

  Scaling parameter for the Laplacian diffusion method (default = 0.7).

- t:

  Time parameter for the Heat Kernel diffusion (default = 1).

- restart_prob:

  Restart probability for the Random Walk with Restart (default = 0.3).

- k:

  Number of seed nodes to select (default = 5).

- normalize:

  Logical; whether to normalize the cumulative diffusion score over the
  target nodes (default = TRUE).

- plot:

  Logical; whether to plot the progression of target diffusion score as
  seeds are selected (default = TRUE).

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

## Value

A list with the following elements:

- selected:

  Character vector of selected seed node names.

- final_target_score:

  Final cumulative (or normalized) diffusion score over the target
  nodes.

- scores_at_each_step:

  Vector of target scores at each greedy selection step.

## Details

This function uses a stepwise greedy heuristic. At each step, it
evaluates the marginal gain in target score from adding each candidate
node to the current seed set, and selects the one with the highest gain.
This is repeated for `k` steps.

Note: The score of each candidate at every step is recomputed via a full
diffusion run, making the function computationally intensive for large
graphs or large `k`.

The target score can be interpreted as either the total or average
diffusion signal received by the target nodes. Normalization helps scale
results across networks of different sizes.

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

## See also

[`network_diffusion`](https://agallinat.github.io/netkit/reference/network_diffusion.md),
[`network_diffusion_with_pvalues`](https://agallinat.github.io/netkit/reference/network_diffusion_with_pvalues.md)

## Examples

``` r
g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
target <- c("1", "2", "3")

res <- greedy_seed_selection(g, target_nodes = target, k = 3, plot = FALSE)
res$selected_seeds
#> [1] "23" "55" "16"
res$scores_at_each_step
#> [1] 0.06180976 0.12123268 0.17697223
```
