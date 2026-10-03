# Assign Vertex and Edge Attributes to an igraph Graph

Adds or updates vertex and edge attributes in an `igraph` object using
user-provided metadata tables. Vertex attributes are matched by the
first column of `nodes_table`, and edge attributes are matched using the
first two columns of `edge_table`, taking graph direction into account.
Only matching nodes and edges are updated. Warnings are issued when
there are unmatched entries.

## Usage

``` r
assign_attributes(
  graph,
  nodes_table = NULL,
  edge_table = NULL,
  overwrite = TRUE
)
```

## Arguments

- graph:

  An `igraph` object or a data frame containing a symbolic edge list.

- nodes_table:

  Optional. A `data.frame` whose first column corresponds to vertex
  names.

- edge_table:

  Optional. A `data.frame` whose first two columns correspond to source
  and target vertices.

- overwrite:

  Logical. If `TRUE`, existing attributes are overwritten. If `FALSE`,
  existing attributes are preserved. Default is `TRUE`.

## Value

An `igraph` object with added or updated attributes.

## Examples

``` r
g <- igraph::make_ring(5)
igraph::V(g)$name <- letters[1:5]

nodes <- data.frame(node = letters[1:5], group = c("a", "a", "b", "b", "b"))
edges <- data.frame(from = c("a", "b"), to = c("b", "c"), weight = c(10, 20))

g <- assign_attributes(g, nodes_table = nodes, edge_table = edges)
igraph::vertex_attr(g, "group")
#> [1] "a" "a" "b" "b" "b"
igraph::edge_attr(g, "weight")
#> [1] 10 20 NA NA NA
```
