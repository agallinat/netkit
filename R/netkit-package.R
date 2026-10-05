#' @keywords internal
#' @aliases netkit-package
#'
#' @section Getting started:
#' Every function accepts either an \pkg{igraph} object or a data.frame edge
#' list, so there is no import step to learn:
#'
#' ```r
#' edges <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"))
#' summarize_graph_metrics(edges)
#' ```
#'
#' Analysis functions return a list with a shared vocabulary -- `result` (a
#' table of per-node or per-step values), `graph` (the input graph with the new
#' values attached as vertex attributes), `method` (what was actually computed,
#' including the thresholds used) and, where there is one, `plot`. Because the
#' annotated graph comes back, the functions chain:
#'
#' ```r
#' hubs <- find_hubs(g, plot = FALSE)
#' mods <- find_modules(hubs$graph, plot = FALSE)   # keeps `is_hub`
#' ```
#'
#' @section Function map:
#' \describe{
#'   \item{Input and annotation}{[assign_attributes()] attaches node and edge
#'     metadata from a data.frame. [netkit-weights] documents how edge weights
#'     are declared and used -- read it before analyzing a weighted graph.}
#'   \item{Topology}{[summarize_graph_metrics()] for one row of global metrics,
#'     [node_metrics()] for the per-node counterpart, [plot_CCDF()] for the
#'     degree distribution and [compare_networks()] for two graphs side by
#'     side.}
#'   \item{Statistical testing}{[null_model()] builds a matched random
#'     ensemble; [metric_significance()] tests observed metrics against it and
#'     [small_worldness()] reports the sigma coefficient.}
#'   \item{Node classification}{[find_hubs()], [find_bottlenecks()] and
#'     [calculate_roles()] for the Guimera-Amaral roles.}
#'   \item{Community structure}{[find_modules()].}
#'   \item{Subnetworks}{[extract_subnetwork()] builds an interpretable
#'     neighborhood around a set of nodes.}
#'   \item{Diffusion}{[network_diffusion()] propagates a signal from seed
#'     nodes, [network_diffusion_with_pvalues()] adds a permutation null,
#'     [prepare_diffusion()] precomputes the kernel for repeated calls and
#'     [greedy_seed_selection()] solves the inverse problem.}
#'   \item{Robustness}{[robustness_analysis()].}
#'   \item{Visualization}{[plot_Net()] is the central renderer;
#'     [highlight_nodes()] and [layout_horizontal_tree()] support it.}
#' }
#'
#' @seealso
#' `vignette("introduction", package = "netkit")` for a worked tour, and
#' `vignette("weighted-networks", package = "netkit")` for what a weight means
#' to each function.
"_PACKAGE"
