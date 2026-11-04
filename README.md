# About `manifoldtest` <img src="man/figures/logo.png" align="right" height="139" alt="" />

`manifoldtest` is an R package to perform independence tests on manifold-valued data. It currently supports data on the $d$-dimensional flat torus, the $d$-dimensional hyperbolic space and the Grassmann manifold $\mathrm{Gr}(n, k)$ of $k$-dimensional linear subspaces of $\mathbb{R}^n$. Tests can also be performed on any other manifold by providing its logarithmic map and geodesic distance, together with a few helper functions (see the user guide).

## Installation

To install the R package from GitHub, run the following command in R:

```r
# install.packages("devtools")
devtools::install_github("gonzalez-delgado/manifoldtest")
```

## User guide

A comprehensive user guide for $\texttt{manifoldtest}$, including tutorials and illustrative examples, is available at <https://gonzalez-delgado.github.io/manifoldtest/>.

## Examples

To reproduce the real-data analyses presented in the paper, see the [`examples/`](examples/) section.