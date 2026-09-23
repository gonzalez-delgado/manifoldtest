# manifoldtest: Independence Testing for Manifold-Valued Data

`manifoldtest` is an R package for testing statistical independence between two samples of data lying on a Riemannian manifold. Each sample is mapped to a Euclidean tangent space via the logarithmic (inverse exponential) map at a chosen base point, optionally transformed to multidimensional ranks bounded in the unit ball via the center-outward distribution function of Hallin et al. (2021), and then tested for independence with the distance covariance test of Székely, Rizzo and Bakirov (2007).

The package currently supports data on the flat torus (e.g. angular or periodic measurements such as dihedral angles), and is designed so that new manifolds can be added from their logarithmic map and distance function.

## Installation

The user can install `manifoldtest` by running the following command:

```r
# install.packages("devtools")
devtools::install_github("gonzalez-delgado/manifoldtest")
```

## Quick example

```r
library(manifoldtest)

# Two samples on the torus, represented in [0,1)^d
X <- matrix(runif(300), nrow = 100, ncol = 3)
Y <- matrix(runif(200), nrow = 100, ncol = 2)

result <- manifold.dcov.test(X, Y, manifold = "torus", R = 500)
result$dcov_test$p.value
```

## On this website

<div class="guide-cards">

  <a class="guide-card" href="reference/">
  <h3>🧩 API documentation</h3>
  <p>Function reference and package documentation.</p>
  </a>

</div>


<div class="citation-box">

<div class="citation-title">Cite this work</div>

If you use `manifoldtest`, please cite:

### The R package

González-Delgado, J. (2026).  
*manifoldtest: Independence testing for manifold-valued data*.  
R package version 0.1.0. <a href="https://github.com/gonzalez-delgado/manifoldtest" target="_blank">GitHub repository</a>

<details class="citation-details">
<summary>BibTeX</summary>

```bibtex
@Manual{manifoldtest2026,
  title   = {manifoldtest: Independence testing for manifold-valued data},
  author  = {González-Delgado, Javier},
  year    = {2026},
  note    = {R package version 0.1.0},
  url     = {https://github.com/gonzalez-delgado/manifoldtest}
}
```
</details> </div>
