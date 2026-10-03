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
  precompute = NULL
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
  this function repeatedly (e.g., in greedy optimization).

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
(optionally normalized) adjacency matrix. For efficiency in iterative
applications, precompute Laplacian and Cholesky decomposition using
[`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md).

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
#>  1 1     0.700 
#>  2 4     0.616 
#>  3 2     0.616 
#>  4 3     0.536 
#>  5 5     0.512 
#>  6 74    0.112 
#>  7 65    0.109 
#>  8 45    0.100 
#>  9 77    0.0926
#> 10 35    0.0902
#> # ℹ 70 more rows

# Reuse a precomputed kernel across repeated calls.
kernel <- prepare_diffusion(g, method = "rwr")
network_diffusion(g, seed_nodes, method = "rwr", precompute = kernel)
#> # A tibble: 80 × 2
#>    node  score
#>    <chr> <dbl>
#>  1 1     0.422
#>  2 3     0.386
#>  3 5     0.383
#>  4 4     0.381
#>  5 2     0.360
#>  6 65    0.161
#>  7 61    0.159
#>  8 45    0.147
#>  9 51    0.137
#> 10 74    0.134
#> # ℹ 70 more rows
```
