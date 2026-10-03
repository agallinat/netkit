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
  normalize = TRUE
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

## Value

A matrix representing the diffusion kernel.

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
#>  1 2     0.707 
#>  2 1     0.663 
#>  3 38    0.0992
#>  4 30    0.0991
#>  5 46    0.0884
#>  6 51    0.0737
#>  7 17    0.0555
#>  8 7     0.0230
#>  9 55    0.0169
#> 10 54    0.0164
#> # ℹ 50 more rows
```
