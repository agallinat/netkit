# Printing netkit objects

netkit's analysis functions return a list holding a `result` table, the
annotated `graph`, a `method` string describing what was computed and,
where there is one, a `plot`. Printing such a list with R's default
method dumps the whole graph and every row of the table, which for a
100-graph
[`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
ensemble runs to well over a thousand lines. These methods print a
summary instead.

## Usage

``` r
# S3 method for class 'netkit_result'
print(x, n = 5, ...)

# S3 method for class 'netkit_null'
print(x, ...)

# S3 method for class 'netkit_null'
x[i]

# S3 method for class 'netkit_kernel'
print(x, ...)
```

## Arguments

- x:

  A netkit object: an analysis result, a
  [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
  ensemble or a
  [`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md)
  kernel.

- n:

  Number of rows of `x$result` to preview. Default 5.

- ...:

  Ignored, present for consistency with the generic.

- i:

  Indices of the null graphs to keep.

## Value

`x`, invisibly. Called for the side effect of printing.

## Details

Only the printing changes. The objects are still plain lists:
`x$result`, `x$graph`, `x$plot` and `x[[i]]` behave exactly as they did,
and `unclass(x)` recovers the default printing.

## Examples

``` r
g <- igraph::sample_gnp(40, 0.1, directed = FALSE)
igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))

find_hubs(g, plot = FALSE)
#> <netkit result: find_hubs()>
#> Method: Hub nodes identified by method: zscore with Degree metric threshold =
#>   3 and Betweenness metric threshold = 1 (unweighted)
#> 
#> $result
#> # A tibble: 40 × 7
#>   node  degree strength betweenness degree_metric betweenness_metric is_hub
#>   <chr>  <dbl>    <dbl>       <dbl>         <dbl>              <dbl> <lgl> 
#> 1 1          3        3      0.0167        -0.448            -0.625  FALSE 
#> 2 2          1        1      0             -2.19             -1.04   FALSE 
#> 3 3          5        5      0.0435         0.572             0.0197 FALSE 
#> 4 4          8        8      0.0746         1.59              0.746  FALSE 
#> 5 5          5        5      0.0662         0.572             0.551  FALSE 
#> # ℹ 35 more rows
#> 
#> $plot   NULL
#> $graph  <igraph> 40 nodes, 83 edges, undirected | vertex attrs: name, is_hub

# The elements are unchanged -- only the printing is summarized.
hubs <- find_hubs(g, plot = FALSE)
head(hubs$result, 3)
#> # A tibble: 3 × 7
#>   node  degree strength betweenness degree_metric betweenness_metric is_hub
#>   <chr>  <dbl>    <dbl>       <dbl>         <dbl>              <dbl> <lgl> 
#> 1 1          3        3      0.0167        -0.448            -0.625  FALSE 
#> 2 2          1        1      0             -2.19             -1.04   FALSE 
#> 3 3          5        5      0.0435         0.572             0.0197 FALSE 

# n = 20 rather than the default 100 null graphs, to keep the example quick.
null_model(g, n = 20, seed = 1)
#> <netkit null ensemble: 20 graphs, model 'rewire'>
#> Each graph has 40 nodes and 83 edges.
#> Use x[[i]] for one graph, or pass the ensemble to metric_significance(null = x).

prepare_diffusion(g, method = "rwr")
#> <netkit diffusion kernel: method 'rwr'>
#> 40 x 40 Laplacian, unweighted edges.
#> Pass to network_diffusion(precompute = x) to skip rebuilding it.
```
