# Edge weights in netkit

How netkit interprets edge weights, and why it asks you to say which
kind you have.

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

## Why the default is to ignore weights

`weights = NULL` ignores edge weights entirely, even when the graph
carries a `weight` attribute. That is deliberate: it makes the
interpretation an explicit choice rather than an accident of which
igraph function happens to be called underneath. When a graph does carry
a `weight` attribute and you pass nothing, netkit warns that it is being
ignored, so the behavior is never silent.

Note that this differs from plain igraph, where the attribute is picked
up automatically.

## Constraints

Weights must be non-negative, finite and non-missing, and must not be
zero for every edge. A zero `"strength"` is allowed and means "no
connection": it maps to an infinite distance, which the path-based
measures treat as unreachable.

Signed attributes – activation versus inhibition, positive versus
negative correlation – are **not** edge weights and are rejected. Take
the absolute value if magnitude is what matters, or keep the sign as a
separate edge attribute for annotation and visualization.
