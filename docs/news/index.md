# Changelog

## netkit 0.0.1.9000 (development version)

Not yet released. The fourth version component (`.9000`) marks this as a
development build; drop it, and update this heading to the release
version, at submission time.

Note the heading must keep a parseable version number. R’s NEWS.md
parser needs one, and a bare `# netkit (development version)` heading
makes `R CMD check --as-cran` report
`Problems with news in 'NEWS.md': No news entries found.`

netkit provides a toolkit for analyzing and visualizing networks, built
on [igraph](https://igraph.org) and returning ggplot2 objects rather
than drawing them, so that diagnostics can be composed and customized by
the caller.

### Features

- **Topology** —
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md)
  produces a one-row table of global metrics (density, diameter, average
  path length, transitivity, assortativity, algebraic connectivity,
  degree entropy, Gini coefficient of degree, modularity).
  [`plot_CCDF()`](https://agallinat.github.io/netkit/reference/plot_CCDF.md)
  draws the complementary cumulative degree distribution on log-log axes
  with optional power-law reference lines, and
  [`compare_networks()`](https://agallinat.github.io/netkit/reference/compare_networks.md)
  puts two graphs side by side with a Kolmogorov-Smirnov test on their
  degree distributions.

- **Node classification** —
  [`find_hubs()`](https://agallinat.github.io/netkit/reference/find_hubs.md)
  (high degree and high betweenness) and
  [`find_bottlenecks()`](https://agallinat.github.io/netkit/reference/find_bottlenecks.md)
  (low degree, high betweenness) support z-score and quantile
  thresholding.
  [`calculate_roles()`](https://agallinat.github.io/netkit/reference/calculate_roles.md)
  implements the Guimerà and Amaral

  2005. R1–R7 connectivity cartography from within-module degree
        z-scores and participation coefficients.

- **Community structure** —
  [`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md)
  wraps nine igraph community detection algorithms with module-size
  filtering and optional subgraph extraction.

- **Diffusion** —
  [`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md)
  propagates influence from seed nodes by Laplacian smoothing, heat
  kernel, or random walk with restart.
  [`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md)
  builds a reusable kernel for repeated runs,
  [`network_diffusion_with_pvalues()`](https://agallinat.github.io/netkit/reference/network_diffusion_with_pvalues.md)
  adds a permutation null parallelized with future.apply, and
  [`greedy_seed_selection()`](https://agallinat.github.io/netkit/reference/greedy_seed_selection.md)
  searches for the seed set that maximizes diffusion onto a target set.

- **Robustness** —
  [`robustness_analysis()`](https://agallinat.github.io/netkit/reference/robustness_analysis.md)
  simulates progressive node removal at random, by degree, by
  betweenness, or by any numeric vertex attribute, and reports
  per-metric area under the curve.

- **Visualization** —
  [`plot_Net()`](https://agallinat.github.io/netkit/reference/plot_Net.md)
  is the central renderer, with
  [`highlight_nodes()`](https://agallinat.github.io/netkit/reference/highlight_nodes.md)
  and
  [`layout_horizontal_tree()`](https://agallinat.github.io/netkit/reference/layout_horizontal_tree.md)
  alongside it.
  [`assign_attributes()`](https://agallinat.github.io/netkit/reference/assign_attributes.md)
  attaches metadata tables to a graph so that results from one function
  can feed the next.

Analysis functions return a named list with a consistent vocabulary —
`plot`, `result`, `graph`, `method` — which is what lets them be
chained.
