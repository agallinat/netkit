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
  verbose = TRUE
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
#>  1 1     0.815       0.0196
#>  2 4     0.746       0.0196
#>  3 3     0.575       0.0196
#>  4 5     0.574       0.0196
#>  5 2     0.541       0.0196
#>  6 22    0.166       0.0196
#>  7 60    0.0949      0.0196
#>  8 31    0.0945      0.0196
#>  9 12    0.0914      0.0196
#> 10 50    0.0876      0.0196
#> # ℹ 50 more rows
```
