# Highlight Nodes in a Network Plot

Highlights selected nodes in an `igraph` object by changing their label,
fill color, and/or outline color. The function modifies node attributes
and visualizes the result using
[`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md).

## Usage

``` r
highlight_nodes(
  graph,
  nodes,
  method = c("label", "fill", "outline"),
  label_color = "darkred",
  highlight_fill_color = "orange",
  highlight_frame_color = "darkblue",
  background_color = "gray",
  ...
)
```

## Arguments

- graph:

  An `igraph` object or a `data.frame` representing a symbolic edge
  list. If a `data.frame`, it should have at least two columns
  specifying source and target nodes.

- nodes:

  A character vector of node names to highlight.

- method:

  Character vector specifying how to highlight the nodes. One or more
  of: `"label"`, `"fill"`, `"outline"`.

- label_color:

  Color for highlighted node labels. Used only if `"label"` is in
  `method`. Default is `"darkred"`.

- highlight_fill_color:

  Fill color for highlighted nodes. Used only if `"fill"` is in
  `method`. Default is `"orange"`.

- highlight_frame_color:

  Outline color for highlighted nodes. Used only if `"outline"` is in
  `method`. Default is `"darkblue"`.

- background_color:

  Fill or frame color for non-highlighted nodes. Default is `"gray"`.

- ...:

  Additional arguments passed to
  [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md).

## Value

Invisibly returns the `igraph` object with updated attributes.

## Examples

``` r
g <- igraph::make_ring(10)
igraph::V(g)$name <- letters[1:10]
highlight_nodes(g, nodes = c("a", "j"), method = "fill")

```
