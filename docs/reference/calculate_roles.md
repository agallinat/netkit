# Calculate Network Roles Based on Within-Module Z-Score and Participation Coefficient

Implements the node role classification system of Guimerà & Amaral
(2005) by calculating the within-module degree z-score and the
participation coefficient for each node in a network. Nodes are assigned
to one of seven role categories (R1–R7) based on their local modular
connectivity.

## Usage

``` r
calculate_roles(
  graph,
  communities = NULL,
  cluster.method = "spinglass",
  plot = TRUE,
  highlight_roles = TRUE,
  hub_z = 2.5,
  label_region = NULL,
  label.size = 12,
  thresholds = NULL,
  weights = NULL,
  weight_type = c("strength", "distance")
)
```

## Arguments

- graph:

  An `igraph` object representing the network, or a data frame
  containing a symbolic edge list in the first two columns. Additional
  columns are considered as edge attributes.

- communities:

  Optional. A community clustering object (as returned by an `igraph`
  clustering function), or a named membership vector. If `NULL`,
  community detection is performed using `cluster.method`.

- cluster.method:

  Character. Clustering algorithm to use if `communities` is `NULL`.
  Default is `"spinglass"`. Passed to
  [`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md).

- plot:

  Logical. Whether to generate a 2D plot of participation
  coefficient (P) vs. within-module z-score (z). Default is `TRUE`.

- highlight_roles:

  Logical. If `TRUE`, the role regions in the z–P plane are shaded for
  visual clarity. Default is `TRUE`.

- hub_z:

  Numeric. Threshold for defining hubs in terms of within-module
  z-score. Default is `2.5`.

- label_region:

  Optional character vector of role labels (e.g., `c("R4", "R7")`)
  indicating which role regions should have their nodes labeled in the
  plot. Default is `NULL`.

- label.size:

  Numeric. Base font size for plot text. Default is `12`.

- thresholds:

  Optional named numeric vector overriding one or more of the
  participation-coefficient boundaries between roles. Names must be
  drawn from `R1_R2`, `R2_R3`, `R3_R4`, `R5_R6` and `R6_R7`; unnamed
  entries, unknown names, values outside `[0, 1]` and non-monotonic sets
  are rejected. `NULL` (default) uses the published values – see
  Details.

- weights:

  Optional edge weights: `NULL` (default) to ignore them, the name of an
  edge attribute, or a numeric vector of length `igraph::ecount(graph)`.
  When supplied, the within-module z-score is computed from strengths
  and the participation coefficient from summed edge strengths per
  module, which is the weighted generalization given in Guimera &
  Amaral's supplementary material. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- weight_type:

  Either `"strength"` (default) or `"distance"`. See
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md).

## Value

A list with five elements:

- `plot`:

  A `ggplot2` object, or `NULL` when `plot = FALSE`. The element is
  always present, so the return shape does not depend on the arguments.

- `result`:

  A data frame with node-level information: node name, module, z-score,
  participation coefficient, and assigned role. It has one row per graph
  vertex. The exception is `cluster.method = "spinglass"`, which can
  only be run on the largest connected component of a disconnected
  graph; vertices outside it have no module and are absent, which raises
  a warning. `z`, `p` and `role` are `NA` for any vertex whose module,
  or whose neighbors' modules, are unknown.

- `graph`:

  The input graph with `module`, `role_z`, `role_p` and `role` attached
  as vertex attributes, so that the classification can be passed
  straight to
  [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md)
  or
  [`robustness_analysis()`](https://agallinat.github.io/netkit/reference/robustness_analysis.md).

- `method`:

  A human-readable description of the community detection used and the
  thresholds actually applied.

- `roles_definitions`:

  A data frame describing the seven role types and their conditions,
  generated from the same thresholds the classifier used.

## Details

If no community structure is provided, modules are automatically
detected using the specified clustering method. The function can
optionally produce a 2D role plot (z vs. P) highlighting the canonical
role regions.

When `communities` is `NULL`, community detection is delegated to
[`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md)
with `min_size = 1`, so no module is discarded for being small and every
vertex receives a role. This matters for correctness as well as
coverage: the participation coefficient of a node is computed from the
module memberships of its neighbors, so dropping a neighbor's module
silently distorts the coefficient of the node that remains.

The node roles are defined as follows, where `hub_z` defaults to 2.5 and
the participation-coefficient boundaries are those published in Guimera
& Amaral (2005):

|     |                                                        |
|-----|--------------------------------------------------------|
| R1  | Ultra-peripheral (non-hub): \\z \< 2.5, P \<= 0.05\\   |
| R2  | Peripheral (non-hub): \\z \< 2.5, 0.05 \< P \<= 0.62\\ |
| R3  | Non-hub connector: \\z \< 2.5, 0.62 \< P \<= 0.80\\    |
| R4  | Non-hub kinless: \\z \< 2.5, P \> 0.80\\               |
| R5  | Provincial hub: \\z \>= 2.5, P \<= 0.30\\              |
| R6  | Connector hub: \\z \>= 2.5, 0.30 \< P \<= 0.75\\       |
| R7  | Kinless hub: \\z \>= 2.5, P \> 0.75\\                  |

Those five numbers are held in one place internally and are used by the
classifier, by the `roles_definitions` table and by the shaded bands of
the diagnostic plot alike, so the three cannot disagree. Override them
with `thresholds` if a different convention is wanted.

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

Guimerà, R., & Amaral, L. A. N. (2005). Functional cartography of
complex metabolic networks. *Nature*, 433(7028), 895–900.
[doi:10.1038/nature03288](https://doi.org/10.1038/nature03288)

## See also

[`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md),
[`igraph::cluster_spinglass()`](https://r.igraph.org/reference/cluster_spinglass.html),
[`igraph::membership()`](https://r.igraph.org/reference/communities.html)

## Examples

``` r
g <- igraph::sample_gnp(80, 0.08, directed = FALSE)
igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))

# "louvain" is used here rather than the "spinglass" default. spinglass cannot
# run on a disconnected graph, so find_modules() falls back to the largest
# connected component and the remaining nodes receive no role.
result <- calculate_roles(g, cluster.method = "louvain", plot = FALSE)
head(result$result)
#> # A tibble: 6 × 5
#>   node  module      z     p role 
#>   <chr>  <int>  <dbl> <dbl> <chr>
#> 1 1          1  1.64  0.617 R2   
#> 2 2          2  0.503 0.5   R2   
#> 3 3          3 -0.504 0.625 R3   
#> 4 4          4 -0.645 0.625 R3   
#> 5 5          2 -1.09  0.625 R3   
#> 6 6          2  0.503 0.5   R2   
result$roles_definitions
#>   Name                Description                      Condition
#> 1   R1 Ultra-peripheral (non-hub)            z < 2.5 & P <= 0.05
#> 2   R2       Peripheral (non-hub) z < 2.5 & 0.05 < P & P <= 0.62
#> 3   R3          Non-hub connector  z < 2.5 & 0.62 < P & P <= 0.8
#> 4   R4            Non-hub kinless              z < 2.5 & P > 0.8
#> 5   R5             Provincial hub            z >= 2.5 & P <= 0.3
#> 6   R6              Connector hub z >= 2.5 & 0.3 < P & P <= 0.75
#> 7   R7                Kinless hub            z >= 2.5 & P > 0.75
```
