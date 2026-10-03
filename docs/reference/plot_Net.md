# Plot an igraph network with customizable node sizes and edge widths

This function plots an `igraph` network graph with options to customize
node size based on vertex degree or a numeric vertex attribute, set node
and edge colors, label nodes, and adjust edge widths (optionally mapped
to edge betweenness centrality).

## Usage

``` r
plot_Net(
  graph,
  label = FALSE,
  color = "#006d77",
  color_ramp = c("blue", "white", "red"),
  NA_color = "gray",
  frame.color = NULL,
  node.size.factor = 1,
  node.degree.map = TRUE,
  edge.color = "#999999",
  edge.width.factor = 1,
  edge.bw.map = TRUE,
  label.color = "#e29578",
  label.size = 1,
  layout = NULL,
  ...
)
```

## Arguments

- graph:

  An `igraph` object representing the network to plot or a data frame
  containing a symbolic edge list in the first two columns. Additional
  columns are considered as edge attributes.

- label:

  Logical, or character vector. If `FALSE`, no node labels are shown. If
  `TRUE`, node labels are set to vertex names or existing labels. If a
  character vector of length equal to the number of vertices, used as
  custom node labels.

- color:

  Character. Color for node fill, or existing vertex attribute name to
  map node fill color. Default is `'#006d77'`.

- color_ramp:

  A character vector of colors (e.g., `c("blue", "white", "red")`) used
  to create a continuous color ramp for mapping numeric node attributes
  to colors. Ignored if a fixed color is provided or if node color is
  categorical. Passed to
  [`colorRampPalette`](https://rdrr.io/r/grDevices/colorRamp.html) to
  interpolate a gradient of colors.

- NA_color:

  Node color for NAs in node color mapping.

- frame.color:

  Character. Color for node frame, or existing vertex attribute name to
  map node frame color. Default is `NULL`.

- node.size.factor:

  Numeric or character. If numeric, acts as a constant scaling factor
  for node size. If character, it should be the name of a numeric vertex
  attribute used to size nodes.

- node.degree.map:

  Logical. If `TRUE` and `node.size.factor` is numeric, node size is
  mapped to vertex degree multiplied by `node.size.factor`. Ignored if
  `node.size.factor` is character.

- edge.color:

  Character. Color of edges. Default is `'#999999'`.

- edge.width.factor:

  Numeric. Factor to scale edge widths. Default is 1.

- edge.bw.map:

  Logical. If `TRUE`, edge widths are mapped to the log-transformed edge
  betweenness centrality. Otherwise, edges have uniform width.

- label.color:

  Character. Color of node labels. Default is `'#e29578'`.

- label.size:

  Numeric. Character expansion factor for label size. Default is 1.

- layout:

  Optional numeric matrix specifying vertex coordinates for layout.

- ...:

  Additional parameters passed to `plot.igraph`.

## Value

Invisibly returns `NULL`. The function produces a plot.

## Examples

``` r
g <- igraph::sample_gnp(50, 0.06, directed = FALSE)
plot_Net(g, label = TRUE, node.size.factor = 2)
#> Vertex attributes 'label' and 'name' are missing, using vertex index for labels.

```
