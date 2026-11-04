#' Hyperbolic space
#'
#' Returns the \code{\link{new_manifold}} object describing the hyperbolic
#' space \eqn{H^d} of constant curvature \eqn{-1}, in the hyperboloid model
#' \deqn{H^d = \{x \in R^{d+1} : -x_0^2 + x_1^2 + \dots + x_d^2 = -1,\ x_0 > 0\}.}
#' Points are represented as rows \eqn{(x_0, x_1, \ldots, x_d)}, so data on
#' \eqn{H^d} have \code{d + 1} columns. This is the object used by
#' \code{\link{manifold.dcov.test}} and \code{\link{choose_base_point}} when
#' \code{manifold = "hyperbolic"}.
#'
#' The exponential map at any point is a global diffeomorphism, so the cut
#' locus is empty and every base point is admissible. Its components are
#' \code{\link{is_in_hyperbolic}}, \code{\link{log_map_hyperbolic}},
#' \code{\link{dist_hyperbolic}}, \code{\link{is_valid_hyperbolic}}, and two
#' base point methods: \code{"auto"} (the default, the closed-form
#' Lorentzian centroid \code{\link{auto_base_point_hyperbolic}}) and
#' \code{"optimize"} (\code{\link{optimize_base_point_hyperbolic}}).
#'
#' @return An object of class \code{"manifold"}.
#'
#' @examples
#' manifold_hyperbolic()
#'
#' @export
manifold_hyperbolic <- function() {

    new_manifold(
        name = "hyperbolic",
        representation = "the hyperboloid model {x in R^(d+1): -x_0^2 + x_1^2 + ... + x_d^2 = -1, x_0 > 0}",
        contains = is_in_hyperbolic,
        log_map = log_map_hyperbolic,
        dist = dist_hyperbolic,
        is_valid = is_valid_hyperbolic,
        base_point = list(

            auto = function(X, ...) {
                list(p = auto_base_point_hyperbolic(X), convergence = NA_integer_)
            },

            optimize = function(X, n_start = 10, maxit = 200, seed = 42, ...) {
                optimize_base_point_hyperbolic(X=X,
                                               n_start=n_start,
                                               maxit=maxit,
                                               seed=seed)
            }
        )
    )
}

# Lorentz boost sending e_0 = (1, 0, ..., 0) to p (inverse = FALSE), or p to
# e_0 (inverse = TRUE). Both are symmetric (d+1) x (d+1) matrices, so they
# act on data stored as rows by right multiplication.
lorentz_boost <- function(p, inverse = FALSE) {

    p0 <- p[1]
    pr <- p[-1]
    s <- if (inverse) -1 else 1
    rbind(c(p0, s * pr),
          cbind(s * pr, diag(length(pr)) + tcrossprod(pr) / (1 + p0)))
}

#' Check that data lies on the hyperboloid
#'
#' Checks whether every row \eqn{x} of \code{X} satisfies
#' \eqn{-x_0^2 + x_1^2 + \dots + x_d^2 = -1} and \eqn{x_0 > 0}, up to a
#' relative tolerance (see \code{\link{manifold_hyperbolic}}).
#'
#' @param X A matrix or data frame of (nominally) hyperbolic data; rows are
#'   samples, with \code{d + 1} columns \eqn{(x_0, \ldots, x_d)}.
#' @param tol Relative tolerance on the hyperboloid equation, scaled by
#'   \eqn{x_0^2}. Defaults to \code{1e-6}.
#'
#' @return Logical: \code{TRUE} if every row of \code{X} lies on the
#'   hyperboloid, \code{FALSE} otherwise.
#'
#' @examples
#' X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
#' is_in_hyperbolic(X)
#'
#' is_in_hyperbolic(matrix(runif(30), nrow = 10, ncol = 3))
#'
#' @export
is_in_hyperbolic <- function(X, tol = 1e-6) {

    X <- as.matrix(X)
    if (ncol(X) < 2 || !is.numeric(X)) return(FALSE)
    x0 <- X[, 1]
    lorentz <- -x0^2 + rowSums(X[, -1, drop = FALSE]^2)
    all(x0 > 0 & abs(lorentz + 1) <= tol * x0^2)
}

#' Map hyperbolic data to the tangent space at a base point
#'
#' Applies the logarithmic map of \eqn{H^d} at base point \code{p}. At
#' \eqn{p = (1, 0, \ldots, 0)} it is
#' \eqn{\log_p(x) = \mathrm{arcsinh}(\|x'\|)\, x' / \|x'\|} with
#' \eqn{x' = (x_1, \ldots, x_d)}, which equals the
#' \eqn{\mathrm{arccosh}(x_0)\, x'/\|x'\|} of the hyperboloid model but is
#' numerically accurate close to \code{p}. For a general \code{p}, the data
#' are first moved by the Lorentz boost sending \code{p} to
#' \eqn{(1, 0, \ldots, 0)}, which identifies the tangent space at \code{p}
#' with \eqn{R^d}. The map is a global diffeomorphism (\eqn{H^d} has no cut
#' locus).
#'
#' @param X A matrix or data frame of points on the hyperboloid; rows are
#'   samples, with \code{d + 1} columns.
#' @param p Numeric vector of length \code{ncol(X)}: the base point on the
#'   hyperboloid.
#'
#' @return A numeric matrix with \code{nrow(X)} rows and \code{d} columns:
#'   the tangent-space coordinates of each row of \code{X}, whose norms are
#'   the geodesic distances to \code{p}.
#'
#' @examples
#' X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
#' log_map_hyperbolic(X, p = c(1, 0, 0))
#'
#' @export
log_map_hyperbolic <- function(X, p) {

    X <- as.matrix(X)
    Y <- X %*% lorentz_boost(p, inverse = TRUE)
    Yr <- Y[, -1, drop = FALSE]
    nr <- sqrt(rowSums(Yr^2))
    Yr * (asinh(nr) / ifelse(nr > 0, nr, 1))
}

#' Map tangent vectors to hyperbolic space
#'
#' Exponential map of \eqn{H^d} at base point \code{p}, the inverse of
#' \code{\link{log_map_hyperbolic}}. Useful to simulate data on \eqn{H^d}.
#'
#' @param V Numeric matrix with \code{d} columns: tangent-space coordinates
#'   at \code{p}, one row per point.
#' @param p Numeric vector of length \code{d + 1}: the base point on the
#'   hyperboloid. Defaults to \eqn{(1, 0, \ldots, 0)}.
#'
#' @return A numeric matrix with \code{nrow(V)} rows and \code{d + 1}
#'   columns: points on the hyperboloid.
#'
#' @examples
#' X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
#' is_in_hyperbolic(X)
#'
#' @export
exp_map_hyperbolic <- function(V, p = c(1, rep(0, ncol(V)))) {

    V <- as.matrix(V)
    r <- sqrt(rowSums(V^2))
    Y <- cbind(cosh(r), V * (sinh(r) / ifelse(r > 0, r, 1)))
    Y %*% lorentz_boost(p)
}

#' Hyperbolic distance from a data set to a point
#'
#' Computes the geodesic distance on \eqn{H^d} between each row of \code{X}
#' and the point \code{p}, i.e. \eqn{\mathrm{arccosh}(-\langle x, p
#' \rangle_L)} with \eqn{\langle x, p \rangle_L = -x_0 p_0 + \sum_i x_i p_i},
#' evaluated as the norm of \code{\link{log_map_hyperbolic}} for numerical
#' accuracy.
#'
#' @param X A matrix or data frame of points on the hyperboloid; rows are
#'   samples, with \code{d + 1} columns.
#' @param p Numeric vector of length \code{ncol(X)}: a point on the
#'   hyperboloid.
#'
#' @return Numeric vector of length \code{nrow(X)}: the hyperbolic distance
#'   from each row of \code{X} to \code{p}.
#'
#' @examples
#' X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
#' dist_hyperbolic(X, c(1, 0, 0))
#'
#' @export
dist_hyperbolic <- function(X, p) {

    sqrt(rowSums(log_map_hyperbolic(X, p)^2))
}

#' Check whether a base point is valid for a hyperbolic data set
#'
#' Checks that \code{p} has one coordinate per column of \code{X} and lies
#' on the hyperboloid. Since \eqn{H^d} has an empty cut locus, any such
#' point is a valid base point.
#'
#' @param X A matrix or data frame of points on the hyperboloid; rows are
#'   samples, with \code{d + 1} columns.
#' @param p Numeric vector of length \code{ncol(X)} giving the base point.
#'
#' @return Logical: \code{TRUE} (invisibly). Stops with an error if
#'   \code{p} has the wrong length or does not lie on the hyperboloid.
#'
#' @examples
#' X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
#' is_valid_hyperbolic(X, p = c(1, 0, 0))
#'
#' @export
is_valid_hyperbolic <- function(X, p) {

    if (length(p) != ncol(X)) {
        stop(sprintf("The length of the base point (%d) mismatches the number of columns of X (%d)!", length(p), ncol(X)))
    }
    if (!is_in_hyperbolic(matrix(p, nrow = 1))) {
        stop("The base point does not lie on the hyperboloid.")
    }

    invisible(TRUE)
}

#' Find a hyperbolic base point via the Lorentzian centroid
#'
#' Closed-form, data-centered base point for \eqn{H^d}: the Euclidean mean
#' \eqn{m} of the rows of \code{X}, rescaled back onto the hyperboloid,
#' \eqn{p = m / \sqrt{-\langle m, m \rangle_L}}. The mean of points on the
#' upper sheet is always future timelike, so \code{p} is well defined. It
#' is the minimizer of the sum of squared Lorentzian distances to the data,
#' and is close to the Fréchet mean when the data are concentrated.
#'
#' @param X A matrix or data frame of points on the hyperboloid; rows are
#'   samples, with \code{d + 1} columns.
#'
#' @return Numeric vector of length \code{ncol(X)}: the chosen base point.
#'
#' @examples
#' X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
#' auto_base_point_hyperbolic(X)
#'
#' @export
auto_base_point_hyperbolic <- function(X) {

    m <- colMeans(as.matrix(X))
    m / sqrt(m[1]^2 - sum(m[-1]^2))
}

#' Find a hyperbolic base point via multi-start optimization
#'
#' Hyperbolic case of \code{\link{optimize_base_point}}: numerically
#' minimizes the total hyperbolic distance (\code{\link{dist_hyperbolic}})
#' from a candidate point \code{p} to every row of \code{X}. The search is
#' unconstrained, over the tangent space at the Lorentzian centroid
#' (\code{\link{auto_base_point_hyperbolic}}), with restarts drawn from a
#' Gaussian scaled to the spread of the data. See
#' \code{\link{choose_base_point}}.
#'
#' @inheritParams optimize_base_point_torus
#' @param X A matrix or data frame of points on the hyperboloid; rows are
#'   samples, with \code{d + 1} columns.
#'
#' @return A list with three elements: \code{p} (the best point found
#'   across all restarts, on the hyperboloid), \code{value} (its total
#'   hyperbolic distance to the rows of \code{X}), and \code{convergence}
#'   (the \code{\link[stats]{optim}} convergence code of the best restart;
#'   \code{0} indicates successful convergence).
#'
#' @examples
#' X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
#' optimize_base_point_hyperbolic(X, n_start = 5)
#'
#' @export
optimize_base_point_hyperbolic <- function(X, n_start = 10, maxit = 200, seed = 42) {

    X <- as.matrix(X)
    d <- ncol(X) - 1
    center <- auto_base_point_hyperbolic(X)
    to_point <- function(v) as.vector(exp_map_hyperbolic(matrix(v, nrow = 1), center))
    spread <- sqrt(mean(dist_hyperbolic(X, center)^2) / d)

    # Optimize over the tangent coordinates at the centroid, which
    # parametrize all of H^d
    optimize_base_point(X = X,
                        dist = function(X, v) dist_hyperbolic(X, to_point(v)),
                        init = function(n) matrix(stats::rnorm(n * d, sd = spread), ncol = d),
                        project = to_point,
                        n_start = n_start,
                        maxit = maxit,
                        seed = seed)
}
