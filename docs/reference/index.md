# Package index

## Graph input and annotation

Getting data into a graph and metadata onto it. Every function in netkit
accepts either an `igraph` object or a data.frame edge list.

- [`assign_attributes()`](https://agallinat.github.io/netkit/reference/assign_attributes.md)
  : Assign Vertex and Edge Attributes to an igraph Graph
- [`netkit-weights`](https://agallinat.github.io/netkit/reference/netkit-weights.md)
  : Edge weights in netkit

## Topology

Global and node-level descriptions of a network’s structure.

- [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md)
  : Summarize Topological Properties of a Graph
- [`node_metrics()`](https://agallinat.github.io/netkit/reference/node_metrics.md)
  : Compute a Table of Node-Level Centrality Metrics
- [`plot_CCDF()`](https://agallinat.github.io/netkit/reference/plot_CCDF.md)
  : Plot Complementary Cumulative Degree Distribution (CCDF)
- [`compare_networks()`](https://agallinat.github.io/netkit/reference/compare_networks.md)
  : Compare Two Networks

## Statistical testing

Whether an observed value is surprising, judged against a matched random
ensemble.

- [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
  : Generate Null-Model Graphs for Significance Testing
- [`metric_significance()`](https://agallinat.github.io/netkit/reference/metric_significance.md)
  : Score an Observed Graph Against a Null Ensemble
- [`small_worldness()`](https://agallinat.github.io/netkit/reference/small_worldness.md)
  : Small-World Coefficients

## Node classification

Picking out structurally distinctive vertices.

- [`find_hubs()`](https://agallinat.github.io/netkit/reference/find_hubs.md)
  : Identify Hub Nodes in a Network
- [`find_bottlenecks()`](https://agallinat.github.io/netkit/reference/find_bottlenecks.md)
  : Identify and Bottleneck Nodes in a Network
- [`calculate_roles()`](https://agallinat.github.io/netkit/reference/calculate_roles.md)
  : Calculate Network Roles Based on Within-Module Z-Score and
  Participation Coefficient

## Community structure

- [`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md)
  : Detect and Visualize Network Modules (Communities)

## Subnetwork extraction

Building an interpretable neighbourhood around a set of nodes.

- [`extract_subnetwork()`](https://agallinat.github.io/netkit/reference/extract_subnetwork.md)
  : Extract an Interpretable Subnetwork Around a Set of Nodes

## Diffusion

Propagating a signal from seed nodes, with a permutation null, and the
inverse problem of choosing seeds to reach a target set.

- [`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md)
  : Perform Network Diffusion from Seed Nodes
- [`network_diffusion_with_pvalues()`](https://agallinat.github.io/netkit/reference/network_diffusion_with_pvalues.md)
  : Perform Network Diffusion from Seed Nodes
- [`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md)
  : Prepare Diffusion Matrix
- [`greedy_seed_selection()`](https://agallinat.github.io/netkit/reference/greedy_seed_selection.md)
  : Greedy Seed Node Selection to Maximize Diffusion Toward Target Nodes

## Robustness

- [`robustness_analysis()`](https://agallinat.github.io/netkit/reference/robustness_analysis.md)
  : Network Robustness Analysis via Node Removal Simulation

## Visualization

[`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md)
is the central renderer; the others delegate to it. These draw with base
graphics rather than returning a ggplot object.

- [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md)
  : Plot an igraph network with customizable node sizes and edge widths
- [`highlight_nodes()`](https://agallinat.github.io/netkit/reference/highlight_nodes.md)
  : Highlight Nodes in a Network Plot
- [`layout_horizontal_tree()`](https://agallinat.github.io/netkit/reference/layout_horizontal_tree.md)
  : Horizontal Tree Layout for Graph Visualization
