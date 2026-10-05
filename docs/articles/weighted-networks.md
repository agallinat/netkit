# Weighted networks: strength or distance?

Most real networks come with a number on every edge: a confidence score
from an interaction database, a correlation from a co-expression matrix,
a count of supporting publications, a reaction rate. This article is
about what netkit does with that number, and why it asks you to say what
kind of number it is.

``` r

library(netkit)
suppressMessages(library(igraph))
```

## The problem

A weight can mean one of two opposite things.

A **strength** means *larger values bind the endpoints more tightly*. A
BioGRID confidence score of 0.9 is a stronger claim than 0.2. A
co-expression correlation of 0.95 is a closer relationship than 0.3.

A **distance** means *larger values push the endpoints further apart*. A
reaction time, a transit cost, a dissimilarity.

Which one you have changes the answer to almost every network question,
in opposite directions. And `igraph` — which netkit is built on — reads
the `weight` edge attribute implicitly, giving it **different meanings
in different functions**: a cost in
[`betweenness()`](https://r.igraph.org/reference/betweenness.html),
[`distances()`](https://r.igraph.org/reference/distances.html),
[`diameter()`](https://r.igraph.org/reference/diameter.html) and
[`mean_distance()`](https://r.igraph.org/reference/distances.html), but
a strength in
[`cluster_louvain()`](https://r.igraph.org/reference/cluster_louvain.html)
and the other community detection algorithms.

Here is what that looks like. The graph is a barbell: two cliques joined
by a single bridge.

``` r

g <- disjoint_union(make_full_graph(4), make_full_graph(4))
V(g)$name <- c(paste0("L", 1:4), paste0("R", 1:4))
g <- add_edges(g, c("L1", "R1"))

# A confidence score: the bridge is a very well-supported interaction, the
# within-clique edges less so.
E(g)$confidence <- 0.1
E(g)$confidence[E(g)["L1" %--% "R1"]] <- 1.0
```

Read as a strength, that bridge is the most reliable edge in the graph,
so crossing it should be *cheap*. Read as a cost, it is the most
expensive edge, so crossing it should *dominate* any path that uses it.

Comparing the two readings’ path lengths directly would be meaningless —
one is measured in units of confidence and the other in units of
1/confidence. What is comparable is the share of a crossing path that
the bridge accounts for:

``` r

bridge_share <- function(type) {
  w <- as_netkit_weights_demo(g, "confidence", type)
  # Cost of getting from one clique to the other, and the bridge's own cost.
  total <- distances(g, v = "L2", to = "R2", weights = w$distance)[1, 1]
  bridge <- w$distance[as.numeric(E(g)["L1" %--% "R1"])]
  c(path_cost = total, bridge_cost = bridge, bridge_share = bridge / total)
}

rbind(strength = bridge_share("strength"),
      distance = bridge_share("distance"))
#>          path_cost bridge_cost bridge_share
#> strength      21.0           1   0.04761905
#> distance       1.2           1   0.83333333
```

Read as a strength, the bridge is the cheapest hop on the route and
accounts for about 5% of it. Read as a cost, the very same number makes
it the most expensive hop, accounting for over 80%. There is no way to
get this right by inspection — both sets of numbers look perfectly
plausible.

## The default is to ignore weights, loudly

netkit’s `weights` argument defaults to `NULL`, which means **ignore
edge weights entirely**, even when the graph carries a `weight`
attribute. When it finds one it is ignoring, it says so:

``` r

E(g)$weight <- E(g)$confidence

res <- tryCatch(
  summarize_graph_metrics(g),
  warning = function(w) conditionMessage(w)
)
res
#> [1] "Graph has an edge attribute 'weight' which is being ignored. Pass weights = \"weight\" to use it, and set weight_type to declare whether it is a 'strength' (larger = more tightly connected) or a 'distance' (larger = further apart)."
```

This differs from plain `igraph`, where the attribute is picked up
automatically. The difference is deliberate: picking it up automatically
is how the three-way inconsistency above stays invisible. Naming the
attribute, or having no `weight` attribute at all, is silent.

## What changes when weights are used

Every function that can use weights takes the same two arguments, and
derives both a strength and a cost vector from them so each metric
receives the one it needs.

``` r

set.seed(1)
net <- sample_gnp(60, 0.08, directed = FALSE)
V(net)$name <- paste0("g", seq_len(vcount(net)))
E(net)$confidence <- runif(ecount(net), 0.05, 1)

unweighted <- summarize_graph_metrics(net)
weighted <- summarize_graph_metrics(net, weights = "confidence")

compare <- rbind(unweighted = unweighted, weighted = weighted)
compare[, c("Avg_degree", "Avg_strength", "Average_path_length",
            "Clustering_coefficient", "Avg_betweenness")]
#>            Avg_degree Avg_strength Average_path_length Clustering_coefficient
#> unweighted   4.366667     4.366667            2.879661             0.05234060
#> weighted     4.366667     2.335531            5.671362             0.04942465
#>            Avg_betweenness
#> unweighted        55.45000
#> weighted          69.41667
```

`Avg_degree` is always the unweighted count and `Avg_strength` the
weighted sum, so the two rows stay comparable. The path-based columns
switch to their weighted definitions, and `Clustering_coefficient`
becomes Barrat’s weighted transitivity.

[`node_metrics()`](https://agallinat.github.io/netkit/reference/node_metrics.md)
makes the same distinction per vertex:

``` r

nm <- node_metrics(net, metrics = c("degree", "strength", "betweenness"),
                   weights = "confidence", normalized = FALSE, plot = FALSE)
head(nm$result)
#> # A tibble: 6 × 4
#>   node  degree strength betweenness
#>   <chr>  <dbl>    <dbl>       <dbl>
#> 1 g1         7     3.22         105
#> 2 g2         5     3.10         142
#> 3 g3         5     2.65          70
#> 4 g4         8     3.45          96
#> 5 g5         2     1.34           6
#> 6 g6         2     1.16          11
```

## Diffusion propagates along strengths

Signal propagation is the one place where there is no ambiguity: a
stronger edge carries more signal. This is also where weights used to be
dropped silently, so it is worth seeing the effect directly.

Back to the barbell, with the bridge almost severed:

``` r

E(g)$w <- 1
g <- set_edge_attr(g, "w", index = E(g)["L1" %--% "R1"], value = 0.001)

far <- paste0("R", 1:4)

open <- network_diffusion(g, "L2", method = "rwr")
#> Warning: Graph has an edge attribute 'weight' which is being ignored. Pass
#> weights = "weight" to use it, and set weight_type to declare whether it is a
#> 'strength' (larger = more tightly connected) or a 'distance' (larger = further
#> apart).
blocked <- network_diffusion(g, "L2", method = "rwr", weights = "w")

c(unweighted = mean(open$score[open$node %in% far]),
  weighted   = mean(blocked$score[blocked$node %in% far]))
#>   unweighted     weighted 
#> 1.617506e-02 3.651383e-05
```

With every edge treated alike, signal from `L2` reaches the far clique
freely. With the bridge at weight 0.001, almost none of it gets through
— which is what a near-severed connection should do.

A strength of exactly zero means “no connection”: it becomes an infinite
distance, which the path-based measures treat as unreachable, and blocks
diffusion completely.

## Seeding from a continuous signal

Related but separate: `seed_weights` lets the diffusion start from
values rather than from set membership. If your seeds came from a
differential-expression result, they are not equally important, and
their magnitudes carry the evidence.

``` r

seeds <- c("g1", "g2", "g3")

binary <- network_diffusion(net, seeds, method = "rwr")
graded <- network_diffusion(net, seeds, method = "rwr",
                            seed_weights = c(g1 = 3.1, g2 = 0.4, g3 = 1.2))

head(data.frame(node = binary$node,
                binary = binary$score,
                graded = graded$score[match(binary$node, graded$node)]))
#>   node    binary    graded
#> 1   g3 0.3930009 0.4814587
#> 2   g1 0.3721936 1.1212384
#> 3   g2 0.3567596 0.1739206
#> 4   g5 0.1741285 0.2463560
#> 5  g11 0.1478400 0.1826112
#> 6  g15 0.1459178 0.4295866
```

## What is not a weight

Signed attributes are rejected rather than coerced:

``` r

E(net)$effect <- sample(c(-1, 1), ecount(net), replace = TRUE)
summarize_graph_metrics(net, weights = "effect")
#> Error: 'weights' contains negative values. Weights must be non-negative; a signed attribute (such as activation/inhibition) is not an edge weight and cannot be used as one.
```

Activation versus inhibition, or a signed correlation, is not an edge
weight — the reciprocal conversion between strength and distance has no
meaning for a negative number. Take the absolute value if magnitude is
what matters, and keep the sign as a separate edge attribute for
annotation and visualization.

## Reusing a diffusion kernel

[`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md)
builds the kernel once so repeated runs can skip it. A kernel is only
valid for the weights it was built from, and netkit checks rather than
trusting you to remember:

``` r

kernel <- prepare_diffusion(net, method = "rwr", weights = "confidence")

# Correct: same weights as the kernel was built with.
head(network_diffusion(net, seeds, method = "rwr", weights = "confidence",
                       precompute = kernel), 3)
#> # A tibble: 3 × 2
#>   node  score
#>   <chr> <dbl>
#> 1 g1    0.393
#> 2 g3    0.378
#> 3 g2    0.362

# Rejected: the kernel describes a different graph than the one requested.
network_diffusion(net, seeds, method = "rwr", precompute = kernel)
#> Error: Pre-computed kernel was built with different edge weights (weighted:strength:131:70.0659407979926:47.5304554513456) than requested (unweighted). Rebuild it with prepare_diffusion(), or pass the same 'weights' and 'weight_type' used to build it.
```

## Summary

| If your edge attribute is | pass |
|----|----|
| a confidence score, correlation, count, interaction score | `weight_type = "strength"` (the default) |
| a cost, a dissimilarity, a transit time | `weight_type = "distance"` |
| a sign (activation/inhibition) | not a weight — use [`abs()`](https://rdrr.io/r/base/MathFun.html), or keep it as annotation |
| not meant to be used at all | nothing; `weights = NULL` is the default |

The one rule worth remembering: **netkit never reads a weight you did
not ask it to read**, and it tells you when a graph carries one you
ignored.
