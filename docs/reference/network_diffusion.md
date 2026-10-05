# Perform Network Diffusion from Seed Nodes

Applies network diffusion techniques to propagate influence from a set
of seed nodes across a graph. Supports Laplacian smoothing, heat
diffusion, and random walk with restart (RWR).

## Usage

``` r
network_diffusion(
  graph,
  seed_nodes,
  method = c("laplacian", "heat", "rwr"),
  alpha = 0.7,
  t = 1,
  restart_prob = 0.3,
  normalize = TRUE,
  precompute = NULL,
  weights = NULL,
  weight_type = c("strength", "distance"),
  seed_weights = NULL
)
```

## Arguments

- graph:

  An `igraph` object representing the network to analyze or a data frame
  containing a symbolic edge list in the first two columns. Additional
  columns are considered as edge attributes. Must have named vertices.

- seed_nodes:

  Character vector of seed node names (must match `V(graph)$name`).

- method:

  Character. Diffusion method to use:

  - `"laplacian"`: Solves the linear system \\(I + \alpha L)^{-1} f_0\\,
    where \\L\\ is the (normalized) graph Laplacian and \\\alpha\\ is a
    smoothing parameter. Internally, a sparse Cholesky decomposition is
    used for efficiency.

  - `"heat"`: Applies the heat diffusion model \\e^{-tL} f_0\\, where
    \\t\\ controls diffusion time. A truncated Taylor expansion is used
    for approximation.

  - `"rwr"`: Random Walk with Restart. Iteratively solves \\f = (1 -
    r)Pf + r f_0\\, where \\P\\ is the transition matrix and \\r\\ is
    the restart probability.

- alpha:

  Damping factor for the Laplacian method. Default is `0.7`.

- t:

  Time parameter for the heat diffusion method. Default is `1`.

- restart_prob:

  Restart probability (usually between 0.3 and 0.7) for the RWR method.
  Default is `0.3`.

- normalize:

  Logical. Whether to normalize the adjacency matrix (symmetric
  normalization for undirected graphs). Default is `TRUE`.

- precompute:

  Optional list of precompute diffusion matrices (e.g., Laplacian,
  Cholesky factor, or transition matrix). Use
  [`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md)
  to generate this object and avoid redundant computations when calling
  this function repeatedly (e.g., in greedy optimization). A kernel
  built with different edge weights than requested is rejected rather
  than reused.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  Signal propagates along edge *strengths*. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- seed_weights:

  Optional numeric vector of initial values for the seeds, replacing the
  default binary indicator. Either named (matched to `seed_nodes` by
  name) or unnamed and parallel to `seed_nodes`. Use this to diffuse
  from a continuous signal – log fold changes, scores, prior
  probabilities – rather than from set membership, which is what the
  propagation literature generally assumes.

## Value

A tibble with two columns, sorted by descending score:

- node:

  Node name

- score:

  Diffusion score representing influence from the seed nodes

## Details

This function allows flexible application of network diffusion
strategies, useful in systems biology (e.g., gene prioritization,
pathway propagation), network analysis, and disease gene discovery. The
underlying matrix operations are based on well-established diffusion
models from graph theory.

For `"rwr"` (random walk with restart), the algorithm iteratively
propagates scores until convergence based on a row-normalized transition
matrix. Recommended method for large networks.

For `"laplacian"` and `"heat"`, the graph Laplacian is computed from the
(optionally normalized) adjacency matrix. For efficiency in iterative
applications, precompute Laplacian and Cholesky decomposition using
[`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md).

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

Köhler S, Bauer S, Horn D, Robinson PN. Walking the interactome for
prioritization of candidate disease genes. *Am J Hum Genet*.
2008;82(4):949–958.
[doi:10.1016/j.ajhg.2008.02.013](https://doi.org/10.1016/j.ajhg.2008.02.013)

Vanunu O, Magger O, Ruppin E, Shlomi T, Sharan R. Associating genes and
protein complexes with disease via network propagation. *PLoS Comput
Biol*. 2010;6(1):e1000641.
[doi:10.1371/journal.pcbi.1000641](https://doi.org/10.1371/journal.pcbi.1000641)

## Examples

``` r
g <- igraph::sample_gnp(80, 0.06, directed = FALSE)
igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
seed_nodes <- igraph::V(g)$name[1:5]

network_diffusion(g, seed_nodes, method = "laplacian")
#> # A tibble: 80 × 2
#>    node   score
#>    <chr>  <dbl>
#>  1 1     0.704 
#>  2 5     0.685 
#>  3 2     0.623 
#>  4 4     0.605 
#>  5 3     0.567 
#>  6 60    0.162 
#>  7 77    0.135 
#>  8 75    0.0839
#>  9 12    0.0812
#> 10 69    0.0805
#> # ℹ 70 more rows

# Reuse a precomputed kernel across repeated calls.
kernel <- prepare_diffusion(g, method = "rwr")
network_diffusion(g, seed_nodes, method = "rwr", precompute = kernel)
#> # A tibble: 80 × 2
#>    node  score
#>    <chr> <dbl>
#>  1 2     0.429
#>  2 3     0.417
#>  3 5     0.398
#>  4 4     0.393
#>  5 1     0.353
#>  6 60    0.300
#>  7 77    0.292
#>  8 12    0.137
#>  9 50    0.134
#> 10 27    0.104
#> # ℹ 70 more rows

# Diffuse along edge strengths rather than treating every edge alike.
igraph::E(g)$confidence <- runif(igraph::ecount(g), 0.1, 1)
network_diffusion(g, seed_nodes, method = "rwr", weights = "confidence")
#> # A tibble: 80 × 2
#>    node  score
#>    <chr> <dbl>
#>  1 2     0.463
#>  2 3     0.422
#>  3 4     0.396
#>  4 5     0.393
#>  5 1     0.344
#>  6 60    0.324
#>  7 77    0.296
#>  8 12    0.154
#>  9 50    0.146
#> 10 75    0.107
#> # ℹ 70 more rows

# Start from a continuous signal instead of set membership.
network_diffusion(g, seed_nodes, method = "rwr",
                  seed_weights = c(2.4, -1.1, 0.7, 3.0, 0.2))
#> # A tibble: 80 × 2
#>    node  score
#>    <chr> <dbl>
#>  1 4     1.05 
#>  2 1     0.824
#>  3 12    0.360
#>  4 3     0.298
#>  5 5     0.243
#>  6 77    0.209
#>  7 36    0.192
#>  8 69    0.177
#>  9 54    0.162
#> 10 38    0.153
#> # ℹ 70 more rows
```
