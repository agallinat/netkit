# Prepare Diffusion Matrix

Prepares and normalizes the diffusion kernel matrix to be used in
network diffusion.

## Usage

``` r
prepare_diffusion(
  graph,
  method = c("laplacian", "heat", "rwr"),
  alpha = 0.7,
  t = 1,
  restart_prob = 0.3,
  normalize = TRUE,
  weights = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph:

  An `igraph` object or a data frame containing a symbolic edge list in
  the first two columns. Additional columns are considered as edge
  attributes.

- method:

  Character string: one of `"laplacian"`, `"heat"`, or `"rwr"`.

- alpha:

  Numeric (used in `"laplacian"`).

- t:

  Time parameter (used in `"heat"`).

- restart_prob:

  Restart probability (used in `"rwr"`).

- normalize:

  Logical. Whether to symmetrically normalize the adjacency matrix.
  Default is `TRUE`.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  Diffusion propagates along edge *strengths*, so a `"distance"` type is
  converted before use. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

## Value

A list with the kernel components: `L` (the Laplacian), `ch` (its
Cholesky factor, for `"laplacian"`), `P` (the transition matrix, for
`"rwr"`), `use_sparse_P`, `method`, and `weights_key` – a digest of the
weights the kernel was built from, which
[`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md)
checks before reusing it.

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
g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))

kernel <- prepare_diffusion(g, method = "laplacian")
kernel$method
#> [1] "laplacian"

# Passing the kernel back in skips rebuilding it on every call.
network_diffusion(g, seed_nodes = c("1", "2"), method = "laplacian",
                  precompute = kernel)
#> # A tibble: 60 × 2
#>    node   score
#>    <chr>  <dbl>
#>  1 2     0.633 
#>  2 1     0.620 
#>  3 57    0.0967
#>  4 49    0.0611
#>  5 44    0.0551
#>  6 35    0.0511
#>  7 13    0.0505
#>  8 39    0.0420
#>  9 48    0.0407
#> 10 11    0.0403
#> # ℹ 50 more rows

# A weighted kernel propagates along edge strengths.
igraph::E(g)$confidence <- runif(igraph::ecount(g), 0.1, 1)
w_kernel <- prepare_diffusion(g, method = "rwr", weights = "confidence")
```
