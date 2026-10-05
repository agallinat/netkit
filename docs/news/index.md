# Changelog

## netkit 0.1.0

First CRAN release.

The “Breaking changes” below are relative to 0.0.1.9000, which was a
development build and was never released, so nothing here changes
behavior anyone could have depended on.

### Breaking changes

- **Edge weights are now ignored unless you ask for them.** Functions
  that can use edge weights take `weights` and `weight_type` arguments,
  and `weights = NULL` (the default) ignores edge weights *even when the
  graph carries a `weight` attribute*, warning that it is doing so. This
  differs from plain , which picks the attribute up automatically — and
  that automatic behavior was the problem: igraph reads a weight as a
  *cost* in
  [`betweenness()`](https://r.igraph.org/reference/betweenness.html),
  [`distances()`](https://r.igraph.org/reference/distances.html),
  [`diameter()`](https://r.igraph.org/reference/diameter.html) and
  [`mean_distance()`](https://r.igraph.org/reference/distances.html),
  but as a *strength* in
  [`cluster_louvain()`](https://r.igraph.org/reference/cluster_louvain.html)
  and the other community algorithms, while netkit’s own matrix code
  ignored it entirely. Attaching a confidence score therefore inverted
  every path-based metric, was read correctly by community detection,
  and vanished from diffusion. See
  [`?"netkit-weights"`](https://agallinat.github.io/netkit/reference/netkit-weights.md).

- **[`network_diffusion_with_pvalues()`](https://agallinat.github.io/netkit/reference/network_diffusion_with_pvalues.md)
  now defaults to a degree-matched permutation null**
  (`null = "degree_matched"`). Diffusion scores are strongly
  degree-dependent, so the previous uniform seed permutation reported
  much of the degree difference between real and permuted seeds as
  significance. This changes every p-value the function produces. Pass
  `null = "uniform"` for the old behavior.

- **[`robustness_analysis()`](https://agallinat.github.io/netkit/reference/robustness_analysis.md)
  no longer re-seeds the random number generator** on its default call.
  It previously ran `set.seed(seed)` unconditionally, so the documented
  default `seed = NULL` became `set.seed(NULL)`, which re-seeds from the
  clock. Results of a default call therefore change, and a seeded script
  is now actually reproducible.
  [`network_diffusion_with_pvalues()`](https://agallinat.github.io/netkit/reference/network_diffusion_with_pvalues.md)
  had the same defect.

- **`compare_networks()$similarity$edge_overlap` is redefined.** It was
  the identical expression to `jaccard_similarity` — one number reported
  twice under two names. It is now the overlap coefficient,
  `|E1 ∩ E2| / min(|E1|, |E2|)`, which is the informative companion to
  Jaccard when two networks differ greatly in size.

- **[`calculate_roles()`](https://agallinat.github.io/netkit/reference/calculate_roles.md)’s
  R2/R3 participation boundary is corrected** from 0.60 to the 0.62
  published by Guimerà and Amaral (2005), and the R5/R6 boundary is now
  0.30 everywhere. Some nodes near those boundaries change role.

- **[`metric_significance()`](https://agallinat.github.io/netkit/reference/metric_significance.md)
  and
  [`small_worldness()`](https://agallinat.github.io/netkit/reference/small_worldness.md)
  take `n` rather than `n_null`**, matching `null_model(n =)`. The old
  name was not just inconsistent: with formals called `null` and
  `n_null` side by side, `metric_significance(g, n = 20)` partially
  matched both and failed with
  `argument 3 matches multiple formal arguments`. With `n` a formal, it
  is an exact match. The `n_null` column of `small_worldness()$result`
  keeps its name.

### New features

- [`node_metrics()`](https://agallinat.github.io/netkit/reference/node_metrics.md)
  — the node-level counterpart to
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md).
  Computes degree, strength, betweenness, closeness, harmonic
  centrality, eigenvector centrality, PageRank, coreness, local
  clustering, Burt’s constraint and eccentricity in one call, and
  attaches each to the graph as a vertex attribute, so
  `robustness_analysis(removal_strategy = "pagerank")` works without
  further steps. The default plot is a Spearman correlation heatmap of
  the metrics, because the usual mistake with such a table is to read
  highly correlated metrics as independent evidence.

- [`extract_subnetwork()`](https://agallinat.github.io/netkit/reference/extract_subnetwork.md)
  — builds an interpretable subnetwork around a set of nodes, by induced
  subgraph, neighborhood expansion (with `max_degree` to exclude
  promiscuous hubs), the union of all shortest paths between seed pairs,
  an approximate minimum Steiner tree (Kou–Markowsky–Berman), or
  diffusion-based expansion reusing the existing kernel.

- [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md),
  [`metric_significance()`](https://agallinat.github.io/netkit/reference/metric_significance.md)
  and
  [`small_worldness()`](https://agallinat.github.io/netkit/reference/small_worldness.md)
  — matched random ensembles and the tests built on them, which turn
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md)
  from descriptive into inferential. Empirical p-values use the
  `(r + 1) / (n + 1)` convention and are never zero.

- Edge weights throughout:
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md),
  [`node_metrics()`](https://agallinat.github.io/netkit/reference/node_metrics.md),
  [`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md),
  [`network_diffusion()`](https://agallinat.github.io/netkit/reference/network_diffusion.md),
  [`network_diffusion_with_pvalues()`](https://agallinat.github.io/netkit/reference/network_diffusion_with_pvalues.md),
  [`greedy_seed_selection()`](https://agallinat.github.io/netkit/reference/greedy_seed_selection.md),
  [`find_hubs()`](https://agallinat.github.io/netkit/reference/find_hubs.md),
  [`find_bottlenecks()`](https://agallinat.github.io/netkit/reference/find_bottlenecks.md),
  [`calculate_roles()`](https://agallinat.github.io/netkit/reference/calculate_roles.md),
  [`find_modules()`](https://agallinat.github.io/netkit/reference/find_modules.md),
  [`robustness_analysis()`](https://agallinat.github.io/netkit/reference/robustness_analysis.md),
  [`compare_networks()`](https://agallinat.github.io/netkit/reference/compare_networks.md),
  [`plot_CCDF()`](https://agallinat.github.io/netkit/reference/plot_CCDF.md)
  and
  [`extract_subnetwork()`](https://agallinat.github.io/netkit/reference/extract_subnetwork.md)
  all accept `weights` and `weight_type`.

- `network_diffusion(seed_weights =)` — diffuse from a continuous signal
  (log fold changes, scores, prior probabilities) rather than from set
  membership.

- `robustness_analysis(removal_strategy = "strength")`.

- [`print()`](https://rdrr.io/r/base/print.html) methods for every
  object netkit returns, so that a result at the console is a summary
  rather than a dump. The default method printed the whole annotated
  graph and every row of the result table; a 100-graph
  [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
  ensemble ran to over 1400 lines. Analysis results now print their
  `method` line, the head of `result`, and a one-line description of
  each remaining element;
  [`null_model()`](https://agallinat.github.io/netkit/reference/null_model.md)
  ensembles and
  [`prepare_diffusion()`](https://agallinat.github.io/netkit/reference/prepare_diffusion.md)
  kernels print their shape. Only printing changes — the objects are
  still plain lists, and
  [`unclass()`](https://rdrr.io/r/base/class.html) restores the old
  output. See
  [`?"netkit-print"`](https://agallinat.github.io/netkit/reference/netkit-print.md).

- A package-level help page:
  [`?netkit`](https://agallinat.github.io/netkit/reference/netkit-package.md)
  gives the function map and the shared return vocabulary.

- `calculate_roles(thresholds =)` — override the
  participation-coefficient boundaries between roles.

### Bug fixes

- [`calculate_roles()`](https://agallinat.github.io/netkit/reference/calculate_roles.md)
  returned `plot`, `roles_definitions` and `result` with no `graph` and
  no `method`, breaking the shared return vocabulary. It now annotates
  the graph with `module`, `role_z`, `role_p` and `role`.

- The five role boundaries were written out three times — in the
  classifier, in the `roles_definitions` table handed to the caller, and
  in the shaded bands of the plot — and had drifted. The R5/R6 boundary
  was 0.30 in the classifier but 0.25 in the other two, so the shaded
  “provincial hub” region disagreed with the classification it
  illustrated. They now come from one place.

- [`network_diffusion_with_pvalues()`](https://agallinat.github.io/netkit/reference/network_diffusion_with_pvalues.md)
  never called the shared input validator, so the data.frame edge list
  it documents failed with igraph’s own error. It was also counting
  seeds before intersecting them with the graph’s vertices, so any
  absent seed made every permuted set larger than the real one,
  inflating the null and biasing every p-value.

- [`compare_networks()`](https://agallinat.github.io/netkit/reference/compare_networks.md)
  crashed on an edgeless graph.
  [`compute_ccdf()`](https://agallinat.github.io/netkit/reference/compute_ccdf.md)
  now returns an empty table for degenerate input, matching
  [`summarize_graph_metrics()`](https://agallinat.github.io/netkit/reference/summarize_graph_metrics.md),
  which reports `NaN` rather than refusing.

- Random walk with restart used an absolute convergence threshold, so
  the precision achieved depended on the magnitude of the seed vector —
  which now matters, since `seed_weights` lets callers set it. The
  threshold is relative to `sum(abs(f0))`, and the previously unbounded
  iteration has a cap.

- [`compare_networks()`](https://agallinat.github.io/netkit/reference/compare_networks.md)
  used [`cat()`](https://rdrr.io/r/base/cat.html) where
  [`message()`](https://rdrr.io/r/base/message.html) was wanted.

### Other

- `find_modules()$module_table` and `robustness_analysis()$summary` gain
  `result` aliases, so the `plot`/`result`/`graph`/`method` return
  vocabulary has no exceptions. The old names are kept and are
  deprecated.

- The test suite has grown from 425 to 894 expectations across 19 files.

## netkit 0.0.1.9000

First development build. Never released.

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
