if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "node", "degree_metric", "betweenness_metric", "module", "step",
    "score", "ccdf", "removed_frac", "lcc_size", "lcc_size_mean",
    "efficiency", "efficiency_mean", "n_components", "n_components_mean",
    ".x", "ccfd", "graph", "ccdf_in", "ccdf_out", "colorRampPalette",
    "role", "target_score", "target_score",
    # calculate_roles(): the z-P role plane. `z` reads as a global only since the
    # per-row classification loop that happened to bind it locally was replaced
    # by vectorized arithmetic; `p` is still bound, by the plot object itself.
    "z",
    # node_metrics(): the correlation heatmap and the ranking facets
    "metric", "metric_x", "metric_y", "correlation", "value", "label",
    # metric_significance(): the faceted null-distribution panel
    "observed"
  ))
}
