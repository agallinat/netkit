# Perform Network Diffusion from Seed Nodes

Applies network diffusion techniques to propagate influence from a set
of seed nodes across a graph. Supports Laplacian smoothing, heat
diffusion, and random walk with restart (RWR).

## Usage

``` r
network_diffusion_with_pvalues(
  graph,
  seed_nodes,
  method = c("laplacian", "heat", "rwr"),
  alpha = 0.7,
  t = 1,
  restart_prob = 0.3,
  normalize = TRUE,
  n_permutations = 1000,
  seed = NULL,
  verbose = TRUE,
  weights = NULL,
  weight_type = c("strength", "distance"),
  seed_weights = NULL,
  null = c("degree_matched", "uniform"),
  match_pool = 10
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
    \\t\\ controls diffusion time. If available, a sparse approximation
    method is used to avoid dense matrix exponential.

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

- n_permutations:

  Integer. Number of permutations to run for empirical p-value
  estimation (default `1000`).

- seed:

  Optional integer for reproducible random number generation. If `NULL`
  (default), seed is not set.

- verbose:

  Logical. If `TRUE` (default), displays a progress bar during
  permutations.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- seed_weights:

  Optional numeric vector of initial seed values, passed to
  [`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md),
  replacing the default binary indicator. Note that the permutation null
  re-uses these magnitudes on the permuted seed sets, so the null tests
  the *position* of the seeds rather than their values.

- null:

  Which permutation null to compare against:

  `"degree_matched"`

  :   The default. Each real seed is replaced by a randomly chosen
      vertex of similar degree, so the permuted sets have the same
      degree profile as the real one.

  `"uniform"`

  :   Permuted seed sets are drawn uniformly from the non-seed vertices,
      matching only the seed count.

  See Details for why the default is degree-matched.

- match_pool:

  Integer. For `null = "degree_matched"`, how many nearest-degree
  candidates each real seed may be replaced by. Default is `10`. Smaller
  values match degree more tightly but give the null less variance;
  `Inf` makes every candidate eligible and so reduces exactly to
  `null = "uniform"`.

## Value

A data frame with two columns:

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
(optionally normalized) adjacency matrix.

## Parallel execution

The permutation null is evaluated with
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html),
which runs under whichever future plan is currently active. This
function does not set a plan itself, so by default permutations run
**sequentially**. To parallelize, set a plan once in your own session
before calling:

    future::plan("multisession", workers = 4)
    res <- network_diffusion_with_pvalues(g, seed_nodes, n_permutations = 1000)
    future::plan("sequential")   # release the workers when finished

Permutations are reproducible regardless of the plan: `seed` is passed
to `future_lapply(future.seed = )`, which generates parallel-safe RNG
streams.

## Choice of null

Diffusion scores are strongly degree-dependent: a high-degree vertex
receives signal from more directions, so it scores highly almost
regardless of where the seeds are. A null that draws seed sets uniformly
therefore compares a real seed set – which in most applications is
enriched for well-studied, highly connected vertices – against a null
composed mostly of poorly connected ones, and reports much of the degree
difference as significance.

`null = "degree_matched"` removes that confound by replacing each real
seed with a randomly chosen vertex drawn from the `match_pool`
candidates closest to it in degree. The resulting p-values answer "are
these seeds unusually well *placed*?" rather than "are these seeds
unusually well *connected*?", which is almost always the intended
question.

How closely the match can be made is limited by the graph. If the seed
set is extreme – the six highest-degree vertices of a scale-free
network, say – there are simply no other vertices of comparable degree
to permute with, and no permutation null can fully correct for degree.
The function warns when the closest available candidates are still far
from the seeds, so the residual confound is visible rather than assumed
away.

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
g <- igraph::sample_gnp(60, 0.08, directed = FALSE)
igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
seed_nodes <- igraph::V(g)$name[1:5]

# n_permutations is reduced from its default of 1000 to keep the example fast;
# use the default or higher for real analyses.
network_diffusion_with_pvalues(g, seed_nodes, method = "laplacian",
                               n_permutations = 50, seed = 1, verbose = FALSE)
#> Warning: package ‘future’ was built under R version 4.4.3
#> # A tibble: 60 × 3
#>    node   score p_empirical
#>    <chr>  <dbl>       <dbl>
#>  1 5     0.770       0.0196
#>  2 4     0.653       0.0196
#>  3 1     0.622       0.0196
#>  4 2     0.617       0.0196
#>  5 3     0.545       0.0196
#>  6 34    0.130       0.0196
#>  7 23    0.0884      0.0196
#>  8 51    0.0873      0.0196
#>  9 59    0.0846      0.0196
#> 10 24    0.0751      0.0196
#> # ℹ 50 more rows
```
