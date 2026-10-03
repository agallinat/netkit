# netkit 0.0.1

First release.

netkit provides a toolkit for analyzing and visualizing networks, built on
[igraph](https://igraph.org) and returning ggplot2 objects rather than drawing
them, so that diagnostics can be composed and customized by the caller.

## Features

* **Topology** — `summarize_graph_metrics()` produces a one-row table of global
  metrics (density, diameter, average path length, transitivity, assortativity,
  algebraic connectivity, degree entropy, Gini coefficient of degree,
  modularity). `plot_CCDF()` draws the complementary cumulative degree
  distribution on log-log axes with optional power-law reference lines, and
  `compare_networks()` puts two graphs side by side with a Kolmogorov-Smirnov
  test on their degree distributions.

* **Node classification** — `find_hubs()` (high degree and high betweenness) and
  `find_bottlenecks()` (low degree, high betweenness) support z-score and
  quantile thresholding. `calculate_roles()` implements the Guimerà and Amaral
  (2005) R1–R7 connectivity cartography from within-module degree z-scores and
  participation coefficients.

* **Community structure** — `find_modules()` wraps nine igraph community
  detection algorithms with module-size filtering and optional subgraph
  extraction.

* **Diffusion** — `network_diffusion()` propagates influence from seed nodes by
  Laplacian smoothing, heat kernel, or random walk with restart.
  `prepare_diffusion()` builds a reusable kernel for repeated runs,
  `network_diffusion_with_pvalues()` adds a permutation null parallelized with
  future.apply, and `greedy_seed_selection()` searches for the seed set that
  maximizes diffusion onto a target set.

* **Robustness** — `robustness_analysis()` simulates progressive node removal at
  random, by degree, by betweenness, or by any numeric vertex attribute, and
  reports per-metric area under the curve.

* **Visualization** — `plot_Net()` is the central renderer, with
  `highlight_nodes()` and `layout_horizontal_tree()` alongside it.
  `assign_attributes()` attaches metadata tables to a graph so that results from
  one function can feed the next.

Analysis functions return a named list with a consistent vocabulary — `plot`,
`result`, `graph`, `method` — which is what lets them be chained.
