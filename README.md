
<!-- README.md is generated from README.Rmd. Please edit that file -->

# netkit

<!-- badges: start -->
<!-- badges: end -->

The goal of **netkit** is to provide a comprehensive and user-friendly
toolkit for network-based data analysis, particularly for biological
systems analysis. It includes tools for flexible graph annotation,
topological analysis, network visualization, and diffusion-based signal
propagation.

Advanced features include a greedy algorithm for reverse diffusion
analysis — predicting optimal seed node sets that maximize signal
propagation to a specified set of target nodes — and tools for
simulating network robustness under targeted or random node removal,
following the framework of [Albert et al.,
2000](https://www.nature.com/articles/35019019). Additionally,
**netkit** implements node role classification based on within-module
and between-module connectivity, as described by [Guimerà & Amaral,
2005](https://www.nature.com/articles/nature03288).

With a special scope to generate high-quality and interpretable figures
suitable for publication, most of the functions generate both tabular
results and diagnostic plots. The package also offers flexible network
visualization options that support node/edge metadata mapping, dynamic
sizing, and layout control.

Altogether, **netkit** is designed to help researchers explore,
interpret, and visualize complex networks with minimal friction and
maximum insight.

## Installation

You can install the development version of **netkit** from GitHub using
the `devtools` or `remotes` package:

``` r
# If not already installed:
install.packages("devtools")

# Install netkit from GitHub
devtools::install_github("agallinat/netkit")
```

## Documentation

Full usage examples and tutorial are available in the package vignettes:

``` r
browseVignettes("netkit")
```

Or online:

- [Introduction to
  netkit](https://agallinat.github.io/netkit/articles/Introduction.html)

## Contributing

Contributions are welcome! If you’d like to report a bug, suggest a
feature, or improve documentation, please open an issue or submit a pull
request at:

<https://github.com/agallinat/netkit/issues>

For larger changes, feel free to open a discussion first.
