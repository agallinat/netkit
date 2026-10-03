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
  label.size = 12
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

## Value

A list with three elements:

- `plot`:

  A `ggplot2` object, or `NULL` when `plot = FALSE`. The element is
  always present, so the return shape does not depend on the arguments.

- `roles_definitions`:

  A data frame describing the seven role types and their conditions.

- `result`:

  A data frame with node-level information: node name, module, z-score,
  participation coefficient, and assigned role. It has one row per graph
  vertex. The exception is `cluster.method = "spinglass"`, which can
  only be run on the largest connected component of a disconnected
  graph; vertices outside it have no module and are absent, which raises
  a warning. `z`, `p` and `role` are `NA` for any vertex whose module,
  or whose neighbours' modules, are unknown.

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
module memberships of its neighbours, so dropping a neighbour's module
silently distorts the coefficient of the node that remains.

The node roles are defined as:

|     |                                                       |
|-----|-------------------------------------------------------|
| R1  | Ultra-peripheral (non-hub): \\z \< 2.5, P \<= 0.05\\  |
| R2  | Peripheral (non-hub): \\z \< 2.5, 0.05 \< P \<= 0.6\\ |
| R3  | Non-hub connector: \\z \< 2.5, 0.6 \< P \<= 0.8\\     |
| R4  | Non-hub kinless: \\z \< 2.5, P \> 0.8\\               |
| R5  | Provincial hub: \\z \>= 2.5, P \<= 0.3\\              |
| R6  | Connector hub: \\z \>= 2.5, 0.3 \< P \<= 0.75\\       |
| R7  | Kinless hub: \\z \>= 2.5, P \> 0.75\\                 |

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
#> 1 1          1  1.64  0.617 R3   
#> 2 2          2  0.503 0.5   R2   
#> 3 3          3 -0.504 0.625 R3   
#> 4 4          4 -0.645 0.625 R3   
#> 5 5          2 -1.09  0.625 R3   
#> 6 6          2  0.503 0.5   R2   
result$roles_definitions
#>   Name                Description                       Condition
#> 1   R1 Ultra-peripheral (non-hub)             z < 2.5 & P <= 0.05
#> 2   R2       Peripheral (non-hub)   z < 2.5 & 0.05 < P & P <= 0.6
#> 3   R3          Non-hub connector    z < 2.5 & 0.6 < P & P <= 0.8
#> 4   R4            Non-hub kinless               z < 2.5 & P > 0.8
#> 5   R5             Provincial hub            z >= 2.5 & P <= 0.25
#> 6   R6              Connector hub z >= 2.5 & 0.25 < P & P <= 0.75
#> 7   R7                Kinless hub             z >= 2.5 & P > 0.75
```
