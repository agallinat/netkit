# netkit 0.1.2

The deferred half of the external 0.1.0 review: the statistical-design,
interpretability and performance items that the 0.1.1 correctness pass left
open. Every behavioural claim below was measured by running the code, and the
performance figures are before-and-after timings on the same machine rather
than estimates.

## New features

* **`compare_networks()` now reports whether the shared edges are more than
  chance gives.** A shared-edge count on its own is not interpretable: two
  large dense networks overlap by arithmetic alone. The new `edge_enrichment`
  element places the observed overlap against the hypergeometric expectation
  over the vertex pairs the two networks could both have reported, giving
  `expected_shared_edges` (\eqn{m_1 m_2 / \binom{n}{2}}), `fold_enrichment` and
  a one-sided `p_value`, plus the counts they are derived from. The universe is
  the *shared* vertex set, and the two draw sizes are each network's distinct
  non-loop edges inside it -- the count and its denominator come from one
  restriction, which is the arithmetic mistake this package has made twice
  before. Reported as its own element rather than as further columns of
  `similarity`, which is documented as a frame of proportions.

* **`greedy_seed_selection()` gains `objective = "coverage"`**, a saturating
  alternative to the default `"total"`. This exists because the default makes
  the greedy search redundant, and that is now stated plainly in the help page:
  diffusion is a linear operator, so the target score of a seed set is the sum
  of its members' individual contributions, the objective is modular, every
  step ranks the candidates identically, and the result is exactly the top `k`
  single-seed contributions. `"coverage"` passes each target's accumulated
  signal through `1 - exp(-x / m_v)` before summing, which is submodular, so
  gains genuinely diminish and the selection spreads across the target set
  instead of concentrating on whichever target is cheapest to reach. Its value
  reads as an effective number of targets reached, bounded by the target count.
  The default is unchanged, so no existing call returns anything different.

* `greedy_seed_selection()` returns `objective_at_each_step` alongside
  `scores_at_each_step`, and a `method` string. `scores_at_each_step` keeps its
  meaning -- the summed diffusion score over the targets -- under both
  objectives, so the two are comparable on one scale.

* **`network_diffusion_with_pvalues()` tests each vertex only against the
  permutations that did not seed it, and reports the null it was tested
  against.** A permutation that happened to seed vertex *v* deposits a large
  score at *v*, for the same reason the real seeds top the real run, and so says
  nothing about whether *v* is unusually well placed with respect to the real
  seed set. Including those permutations inflated the upper tail of each
  vertex's own null and made its p-value conservative -- hardest on the vertices
  drawn most often, so not a uniform shift. **This changes every p-value the
  function produces.**

  The returned table gains `null_mean`, `null_sd`, `z` and `n_valid`, and is now
  ordered by `p_empirical` and then by descending `z`. The z-score is the point:
  the p-value floor is `1 / (n_valid + 1)`, so at the default 1000 permutations
  every strong vertex reports 0.000999 and the ordering among them was
  previously arbitrary. The help page gains an "Interpreting the p-values"
  section stating that the p-values are one-sided, that `n_valid` sets each
  vertex's attainable floor, and that the real seeds sit at the top of the
  ranking by construction and should be dropped before interpreting anything.

  One interaction is worth knowing and is documented under "Choice of null": the
  degree-matched null is no longer uniformly more conservative than the uniform
  one, node count for node count. Its bins restrict which vertices are ever
  seeded, so it draws each of them more often -- up to 10 of 60 permutations
  against 5 for the uniform null -- and the exclusion therefore removes a larger
  share of its nulls' upper tails. Both remain valid; the matched one is
  conditioned more strongly. A test that asserted the old count relationship has
  been replaced by one asserting the thing it was a proxy for, namely that the
  matched null's variance has not collapsed.

* **`greedy_seed_selection()` reports a degree-matched gain z-score.** The new
  `result` table carries one row per selected seed with its `degree`, the
  marginal `gain` that pick bought, and `gain_z` -- that gain standardized
  against the gains of the candidates closest to it in degree. This separates
  the two reasons a vertex can look good: a `gain_z` near zero means the vertex
  is typical of its degree class and was picked *because of* its degree, while a
  large one means it outperforms others as poorly connected as itself. When a
  selection fills up with leaves, this is what says whether any of them is
  individually worth attention. The comparison pool is taken as a whole degree
  radius rather than the *n* nearest, so ties are never resolved by vertex
  index -- the mistake that made the diffusion null's degree-matched sampler
  degenerate in 0.1.1.

* **`metric_significance()` and `small_worldness()` evaluate their null
  ensembles in parallel.** Both recomputed every metric for every null graph
  serially, reported at roughly 19 s per graph at 6000 vertices against a
  default of `n = 100`. They now share one internal summarizer running on
  [future.apply::future_lapply()]. Measured on a 2500-vertex, 15000-edge graph
  with 12 nulls: 9.9 s to 3.3 s on four workers, for identical results.
  Following the convention the diffusion null already set, neither function
  calls `future::plan()` -- the caller sets the plan and the default stays
  sequential.

## Behaviour changes

* **`find_hubs()` and `find_bottlenecks()` warn when they select nothing.**
  Both used to return an empty table in silence, which conflates three
  different situations, and the warning now says which one occurred because the
  remedy differs for each: a metric with no variance, a threshold no vertex can
  reach, or two criteria each satisfied by some vertex but never by the same
  one.

  The middle case is the one reported from the correlation-network analyses, and
  it is not specific to dense graphs: a z-score threshold is only reachable if
  the metric has a tail to reach it. The `"zscore"` defaults of 3 and 1 are
  calibrated for a scale-free degree distribution, and on a near-Gaussian one
  the largest attainable z-score for a few hundred vertices is around 2.5. The
  default stays `"zscore"` deliberately. `"quantile"` was considered and
  rejected: it cannot return nothing, which sounds like an improvement and is
  not -- it selects a fixed *fraction* of the graph whether or not the graph has
  any hub structure, marking every vertex of a ring a hub. See the new
  "Choosing a method" section of `?find_hubs` for the trade-off. (Measuring
  that comparison also turned up the zero-betweenness defect below, which was
  the more serious half of what `"quantile"` was doing.)

* **`find_hubs()` and `find_bottlenecks()` no longer return an all-`NA`
  classification on a regular graph.** When every vertex carries the same
  degree the standard deviation is zero, `scale()` returns `NaN` throughout, and
  every comparison against it was `NA` -- so `is_hub` was `NA` for every vertex
  and was written onto the graph as an all-`NA` vertex attribute, which
  downstream filters silently dropped. The standardized columns stay `NA`, which
  is what an undefined z-score should report, but the classification is now
  `FALSE`: no vertex is unusual when none differs from any other.

* **`find_hubs()` and `find_bottlenecks()` no longer select vertices whose
  betweenness is exactly zero.** Both are defined by *high* betweenness, and
  zero is not high. The guard matters because betweenness is heavily
  zero-inflated on a sparse graph: 87% of the vertices of a 2000-vertex
  scale-free graph have betweenness 0, which puts the 75th percentile *at* zero,
  and the inclusive `>=` comparison then admitted all of them. The default
  `find_bottlenecks(method = "quantile")` call returned **1731 of 2000 vertices**
  as bottlenecks, every one of them lying on no shortest path at all -- the
  definitional opposite of a bottleneck. No choice of quantile defends against
  this, because the quantile of a mostly-zero distribution is zero. The effect
  on `find_hubs()` is nil on every graph family tested, a high-degree vertex
  normally having paths through it; the guard makes that a contract rather than
  a regularity.

* `greedy_seed_selection()` warns when more than half its selection has degree
  1. Under symmetric normalization an edge into a degree-1 vertex carries the
  largest weight any edge to that target can, so a leaf hanging off a target is
  the most efficient injection point in the graph. That is arithmetically
  correct and rarely actionable, and nothing in the returned scores hinted at
  it.

* `greedy_seed_selection()` caps `k` at the number of candidates instead of
  selecting past the end of the candidate set, and errors on an empty candidate
  set.

## Performance

* **The diffusion permutation null is solved as a multiple-right-hand-side
  problem.** Permutations are grouped into batches and each batch is solved in
  one call, so a sparse Cholesky pays its per-call overhead once per batch
  instead of once per permutation. The permuted seed sets are drawn exactly as
  before -- one `future_lapply` element and so one RNG stream per permutation --
  so for a given `seed` the p-values are unchanged; this was verified
  bit-for-bit on all three methods.

  Two further findings came out of profiling rather than from reading the code.
  The diffusion kernel was being built **twice** per call: the real-score
  diffusion ran before the kernel existed and silently built one of its own,
  which for `"laplacian"` meant a second sparse Cholesky factorization of the
  whole graph whose result was discarded. And the batch right-hand side has to
  be a `dgeMatrix`: a base `matrix` takes a markedly slower path through
  Matrix's sparse product, 3.4 s against 0.6 s for the same 30 random-walk
  iterations over 166 columns, for bit-identical results.

  End to end on a 6000-vertex, 78000-edge graph with 200 permutations:
  `"laplacian"` 6.8 s to 3.4 s, `"heat"` 1.4 s to 0.6 s, `"rwr"` 1.4 s to
  0.6 s. The batching saving is proportional to `n_permutations`, so it matters
  at the documented default of 1000 and is negligible at 50; what dominates at
  small permutation counts is the one-time factorization, which is why building
  it once was the larger win there.

* **`greedy_seed_selection()` computes every diffusion once.** Because the
  objective is additive, the entire search reduces to arithmetic on one
  targets-by-candidates matrix of single-seed contributions, built in a single
  batched pass. The previous implementation ran a full diffusion for every
  candidate at every step, recomputing at step 3 what it had already computed
  at step 1. On a 2000-vertex graph with 4 targets and `k = 8` that is 15968
  diffusions and 15.7 s, against 0.37 s now, for the same selection and the
  same scores to within 3e-16.

* The diffusion maths now has a single implementation, `diffuse_columns()`,
  which takes any number of seed vectors at once; `network_diffusion()` calls it
  with one column. Batching is bit-identical for `"laplacian"` and `"heat"`.
  `"rwr"` is the exception: its iteration stops when the slowest column of a
  batch has converged, so the others are carried marginally past their own
  stopping point and end up slightly *closer* to the fixed point. The
  difference is bounded by the convergence tolerance, and is asserted per column
  rather than assumed.

* `compare_networks()` builds its canonical edge keys with vectorized
  `pmin()`/`pmax()` rather than an `apply()` closure per edge, and resolves
  vertex identifiers itself instead of letting `as_edgelist()` choose -- for an
  unnamed graph it returns doubles, and pasting those goes through scientific
  notation above 1e5, which would key one vertex two ways.

# netkit 0.1.1

A correctness release. Every item below was reported from external analysis
runs against 0.1.0 and reproduced here before being fixed; each has a
regression test.

## Bug fixes

* **`network_diffusion_with_pvalues()` matched each node against the wrong
  null distribution.** `network_diffusion()` returns its nodes sorted by
  descending score, so the real scores and every permutation came back in a
  different vertex order. The permutation vectors were combined with
  `do.call(cbind, ...)`, which binds by position and takes its row names from
  the first element only, so each node was tested against the distribution of
  the k-th *ranked* permuted score rather than its own. The symptom was
  bimodal p-values piled at the floor and at 1: on a 200-vertex scale-free
  graph, 135 of 200 nodes sat at the minimum attainable p-value and the rank
  correlation between score and p-value was -0.50, against -0.92 once fixed.
  **This changes every p-value the function has ever produced.**

* **The degree-matched sampler produced an anti-conservative null.** Each seed
  drew from the `match_pool` nearest-degree candidates chosen with `order()`,
  which breaks ties by vertex index, so seeds of similar degree received
  near-identical pools and then drew from them without replacement. Across 200
  permutations of a 5-seed set on a 150-vertex graph the sampler touched 10 of
  145 candidates. A null that narrow has almost no variance, and the
  degree-matched null consequently reported *more* significant nodes than the
  uniform null it exists to be more conservative than. Candidates are now
  grouped into adaptive degree bins, ties are broken at random, and seeds are
  processed in a random order. `match_pool` now defaults to `NULL`, meaning
  `10 * length(seed_nodes)`: the pool has to scale with the seed count, because
  every permutation consumes `length(seed_nodes)` of it. Setting it small
  enough to re-enter the degenerate regime now warns. This defect was only
  visible once the row alignment above was corrected.

* **`null_model()` recognized edge weights only under the literal attribute
  name `weight`.** A graph annotated with any other name — `corr`,
  `confidence`, `combined_score` — produced unweighted null graphs, and
  `metric_significance()` and `small_worldness()` then failed with
  `Edge attribute 'weight' not found` naming an attribute the caller had never
  mentioned. `null_model()` gains a `weights` argument accepting an attribute
  name or a numeric vector, and the two callers forward theirs to it. The
  permuted values are always written onto the nulls as `weight`, as the
  documentation already claimed.

* **Rewired null graphs inherited the observed graph's edge attributes**, now
  bound to edges that had since been moved. Those values were neither the
  observed annotation nor a permutation of one, but a downstream weighted
  metric read them as data. Null graphs now carry no edge attribute except the
  permuted `weight`, which also makes the three models agree on what they
  return.

* **`plot_Net()` resolved `frame.color` by testing `color`.** Mapping `color`
  to a vertex attribute while passing a literal `frame.color` looked up an
  attribute named, for example, `"black"`, got `NULL` and silently drew no
  frame; the mirror case passed the literal attribute name to `plot.igraph` as
  a colour and errored. Neither combination worked.

* **`Modularity` in `summarize_graph_metrics()` was not reproducible.** Louvain
  randomizes its vertex visiting order, so repeated calls on one graph returned
  values spanning roughly 0.356 to 0.366 — larger than many of the differences
  the metric gets used to argue about. It is now computed under a fixed seed,
  applied locally so that the caller's random stream is restored afterwards.
  The new `seed` argument controls it; `seed = NULL` restores the old
  behaviour.

* **`find_modules(method = "leiden")` put every node in its own community on
  weighted graphs.** It used \pkg{igraph}'s default CPM objective at resolution
  1, under which no grouping beats a singleton once the weights fall below 1;
  `min_size` then discarded all of them, leaving an empty result table with
  nothing to indicate why. The objective now defaults to `"modularity"`, and
  `objective_function`, `resolution` and `n_iterations` are exposed. Supplying
  them to another method warns rather than ignoring them.

## New features

* `find_modules()` returns `modularity`, the \eqn{Q} of the partition it found,
  from `igraph::modularity()` and computed before `min_size` filtering. No
  method exposed it before, so comparing two partitions of one graph meant
  recomputing it by hand.

* `find_modules()` warns when `min_size` filtering leaves no modules at all,
  rather than returning a zero-row table that looks the same whether the graph
  has no community structure or the algorithm was misconfigured.

## Documentation

* `greedy_seed_selection()`: the return element is `selected_seeds`, not
  `selected` as documented. `normalize` is passed to the adjacency
  normalization in `network_diffusion()`; it was documented as normalizing the
  cumulative score over the target nodes, which it has never done.

* `network_diffusion_with_pvalues()` gains a "Choice of null" discussion of why
  matching degree more tightly is not better, and how an over-tight
  `match_pool` makes the p-values anti-conservative.

# netkit 0.1.0

First CRAN release.

The "Breaking changes" below are relative to 0.0.1.9000, which was a
development build and was never released, so nothing here changes behavior
anyone could have depended on.

## Breaking changes

* **Edge weights are now ignored unless you ask for them.** Functions that can
  use edge weights take `weights` and `weight_type` arguments, and `weights =
  NULL` (the default) ignores edge weights *even when the graph carries a
  `weight` attribute*, warning that it is doing so. This differs from plain
  \pkg{igraph}, which picks the attribute up automatically — and that automatic
  behavior was the problem: igraph reads a weight as a *cost* in
  `betweenness()`, `distances()`, `diameter()` and `mean_distance()`, but as a
  *strength* in `cluster_louvain()` and the other community algorithms, while
  netkit's own matrix code ignored it entirely. Attaching a confidence score
  therefore inverted every path-based metric, was read correctly by community
  detection, and vanished from diffusion. See `?"netkit-weights"`.

* **`network_diffusion_with_pvalues()` now defaults to a degree-matched
  permutation null** (`null = "degree_matched"`). Diffusion scores are strongly
  degree-dependent, so the previous uniform seed permutation reported much of
  the degree difference between real and permuted seeds as significance. This
  changes every p-value the function produces. Pass `null = "uniform"` for the
  old behavior.

* **`robustness_analysis()` no longer re-seeds the random number generator** on
  its default call. It previously ran `set.seed(seed)` unconditionally, so the
  documented default `seed = NULL` became `set.seed(NULL)`, which re-seeds from
  the clock. Results of a default call therefore change, and a seeded script is
  now actually reproducible. `network_diffusion_with_pvalues()` had the same
  defect.

* **`compare_networks()$similarity$edge_overlap` is redefined.** It was the
  identical expression to `jaccard_similarity` — one number reported twice under
  two names. It is now the overlap coefficient, `|E1 ∩ E2| / min(|E1|, |E2|)`,
  which is the informative companion to Jaccard when two networks differ greatly
  in size.

* **`calculate_roles()`'s R2/R3 participation boundary is corrected** from 0.60
  to the 0.62 published by Guimerà and Amaral (2005), and the R5/R6 boundary is
  now 0.30 everywhere. Some nodes near those boundaries change role.

* **`metric_significance()` and `small_worldness()` take `n` rather than
  `n_null`**, matching `null_model(n =)`. The old name was not just
  inconsistent: with formals called `null` and `n_null` side by side,
  `metric_significance(g, n = 20)` partially matched both and failed with
  `argument 3 matches multiple formal arguments`. With `n` a formal, it is an
  exact match. The `n_null` column of `small_worldness()$result` keeps its
  name.

## New features

* `node_metrics()` — the node-level counterpart to `summarize_graph_metrics()`.
  Computes degree, strength, betweenness, closeness, harmonic centrality,
  eigenvector centrality, PageRank, coreness, local clustering, Burt's
  constraint and eccentricity in one call, and attaches each to the graph as a
  vertex attribute, so `robustness_analysis(removal_strategy = "pagerank")`
  works without further steps. The default plot is a Spearman correlation
  heatmap of the metrics, because the usual mistake with such a table is to read
  highly correlated metrics as independent evidence.

* `extract_subnetwork()` — builds an interpretable subnetwork around a set of
  nodes, by induced subgraph, neighborhood expansion (with `max_degree` to
  exclude promiscuous hubs), the union of all shortest paths between seed pairs,
  an approximate minimum Steiner tree (Kou–Markowsky–Berman), or diffusion-based
  expansion reusing the existing kernel.

* `null_model()`, `metric_significance()` and `small_worldness()` — matched
  random ensembles and the tests built on them, which turn
  `summarize_graph_metrics()` from descriptive into inferential. Empirical
  p-values use the `(r + 1) / (n + 1)` convention and are never zero.

* Edge weights throughout: `summarize_graph_metrics()`, `node_metrics()`,
  `prepare_diffusion()`, `network_diffusion()`,
  `network_diffusion_with_pvalues()`, `greedy_seed_selection()`, `find_hubs()`,
  `find_bottlenecks()`, `calculate_roles()`, `find_modules()`,
  `robustness_analysis()`, `compare_networks()`, `plot_CCDF()` and
  `extract_subnetwork()` all accept `weights` and `weight_type`.

* `network_diffusion(seed_weights =)` — diffuse from a continuous signal (log
  fold changes, scores, prior probabilities) rather than from set membership.

* `robustness_analysis(removal_strategy = "strength")`.

* `print()` methods for every object netkit returns, so that a result at the
  console is a summary rather than a dump. The default method printed the whole
  annotated graph and every row of the result table; a 100-graph `null_model()`
  ensemble ran to over 1400 lines. Analysis results now print their `method`
  line, the head of `result`, and a one-line description of each remaining
  element; `null_model()` ensembles and `prepare_diffusion()` kernels print
  their shape. Only printing changes — the objects are still plain lists, and
  `unclass()` restores the old output. See `?"netkit-print"`.

* A package-level help page: `?netkit` gives the function map and the shared
  return vocabulary.

* `calculate_roles(thresholds =)` — override the participation-coefficient
  boundaries between roles.

## Bug fixes

* `calculate_roles()` returned `plot`, `roles_definitions` and `result` with no
  `graph` and no `method`, breaking the shared return vocabulary. It now
  annotates the graph with `module`, `role_z`, `role_p` and `role`.

* The five role boundaries were written out three times — in the classifier, in
  the `roles_definitions` table handed to the caller, and in the shaded bands of
  the plot — and had drifted. The R5/R6 boundary was 0.30 in the classifier but
  0.25 in the other two, so the shaded "provincial hub" region disagreed with
  the classification it illustrated. They now come from one place.

* `network_diffusion_with_pvalues()` never called the shared input validator, so
  the data.frame edge list it documents failed with igraph's own error. It was
  also counting seeds before intersecting them with the graph's vertices, so any
  absent seed made every permuted set larger than the real one, inflating the
  null and biasing every p-value.

* `compare_networks()` crashed on an edgeless graph. `compute_ccdf()` now
  returns an empty table for degenerate input, matching
  `summarize_graph_metrics()`, which reports `NaN` rather than refusing.

* Random walk with restart used an absolute convergence threshold, so the
  precision achieved depended on the magnitude of the seed vector — which now
  matters, since `seed_weights` lets callers set it. The threshold is relative to
  `sum(abs(f0))`, and the previously unbounded iteration has a cap.

* `compare_networks()` used `cat()` where `message()` was wanted.

## Other

* `find_modules()$module_table` and `robustness_analysis()$summary` gain
  `result` aliases, so the `plot`/`result`/`graph`/`method` return vocabulary has
  no exceptions. The old names are kept and are deprecated.

* The test suite has grown from 425 to 894 expectations across 19 files.

# netkit 0.0.1.9000

First development build. Never released.

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
