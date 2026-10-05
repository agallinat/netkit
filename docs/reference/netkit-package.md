# netkit: Network Analysis and Visualization Toolkit

A collection of tools for analyzing and visualizing networks. Summarizes
global topology and per-node centrality, identifies hubs, bottlenecks
and community structure, and classifies nodes into connectivity roles
following Guimera and Amaral (2005)
[doi:10.1038/nature03288](https://doi.org/10.1038/nature03288) .
Extracts interpretable subnetworks around a set of nodes, including
approximate Steiner trees after Kou, Markowsky and Berman (1981)
[doi:10.1007/BF00288961](https://doi.org/10.1007/BF00288961) .
Propagates influence from seed nodes using Laplacian smoothing,
heat-kernel and random-walk-with-restart diffusion, with degree-matched
permutation nulls, and simulates network robustness under progressive
node removal. Observed metrics can be tested against degree-preserving
random ensembles. Edge weights are supported throughout, with an
explicit declaration of whether a weight means strength or distance.
Matrix operations use sparse representations throughout so that the
methods remain usable on large graphs.

## Getting started

Every function accepts either an igraph object or a data.frame edge
list, so there is no import step to learn:

    edges <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"))
    summarize_graph_metrics(edges)

Analysis functions return a list with a shared vocabulary – `result` (a
table of per-node or per-step values), `graph` (the input graph with the
new values attached as vertex attributes), `method` (what was actually
computed, including the thresholds used) and, where there is one,
`plot`. Because the annotated graph comes back, the functions chain:

    hubs <- find_hubs(g, plot = FALSE)
    mods <- find_modules(hubs$graph, plot = FALSE)   # keeps `is_hub`

## Function map

- Input and annotation:

  [`assign_attributes()`](https://agallinat.github.io/netkit/reference/assign_attributes.md)
  attaches node and edge metadata from a data.frame.
  [netkit-weights](https://agallinat.github.io/netkit/reference/netkit-weights.md)
  documents how edge weights are declared and used – read it before
  analyzing a weighted graph.

- Topology:

  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md)
  for one row of global metrics,
  [`node_metrics()`](https://agallinat.github.io/netkit/reference/node_metrics.md)
  for the per-node counterpart,
  [`plot_CCDF()`](https://agallinat.github.io/netkit/reference/plot_CCDF.md)
  for the degree distribution and
  [`compare_networks()`](https://agallinat.github.io/netkit/reference/compare_networks.md)
  for two graphs side by side.

- Statistical testing:

  [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
  builds a matched random ensemble;
  [`metric_significance()`](https://agallinat.github.io/netkit/reference/metric_significance.md)
  tests observed metrics against it and
  [`small_worldness()`](https://agallinat.github.io/netkit/reference/small_worldness.md)
  reports the sigma coefficient.

- Node classification:

  [`find_hubs()`](https://agallinat.github.io/netkit/reference/find_hubs.md),
  [`find_bottlenecks()`](https://agallinat.github.io/netkit/reference/find_bottlenecks.md)
  and
  [`calculate_roles()`](https://agallinat.github.io/netkit/reference/calculate_roles.md)
  for the Guimera-Amaral roles.

- Community structure:

  [`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md).

- Subnetworks:

  [`extract_subnetwork()`](https://agallinat.github.io/netkit/reference/extract_subnetwork.md)
  builds an interpretable neighborhood around a set of nodes.

- Diffusion:

  [`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md)
  propagates a signal from seed nodes,
  [`network_diffusion_with_pvalues()`](https://agallinat.github.io/netkit/reference/network_diffusion_with_pvalues.md)
  adds a permutation null,
  [`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md)
  precomputes the kernel for repeated calls and
  [`greedy_seed_selection()`](https://agallinat.github.io/netkit/reference/greedy_seed_selection.md)
  solves the inverse problem.

- Robustness:

  [`robustness_analysis()`](https://agallinat.github.io/netkit/reference/robustness_analysis.md).

- Visualization:

  [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md)
  is the central renderer;
  [`highlight_nodes()`](https://agallinat.github.io/netkit/reference/highlight_nodes.md)
  and
  [`layout_horizontal_tree()`](https://agallinat.github.io/netkit/reference/layout_horizontal_tree.md)
  support it.

## See also

[`vignette("introduction", package = "netkit")`](https://agallinat.github.io/netkit/articles/introduction.md)
for a worked tour, and
[`vignette("weighted-networks", package = "netkit")`](https://agallinat.github.io/netkit/articles/weighted-networks.md)
for what a weight means to each function.

## Author

**Maintainer**: Alex Gallinat <alex.gaoc@gmail.com>
