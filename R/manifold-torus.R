#' The flat torus manifold
#'
#' Returns the \code{\link{new_manifold}} object describing the flat torus
#' \eqn{T^d}, with points represented as rows in \eqn{[0,1)^d}. This is the
#' object used by \code{\link{manifold.dcov.test}} and
#' \code{\link{choose_base_point}} when \code{manifold = "torus"}.
#'
#' Its components are \code{\link{is_in_torus}}, \code{\link{log_map_torus}},
#' \code{\link{dist_torus}}, \code{\link{is_valid_torus}}, and two base
#' point methods: \code{"auto"} (the default,
#' \code{\link{auto_base_point_torus}}) and \code{"optimize"}
#' (\code{\link{optimize_base_point_torus}}).
#'
#' @return An object of class \code{"manifold"}.
#'
#' @examples
#' manifold_torus()
#'
#' @export
manifold_torus <- function() {

    new_manifold(
        name = "torus",
        representation = "[0,1)^d",
        contains = is_in_torus,
        log_map = log_map_torus,
        dist = dist_torus,
        is_valid = is_valid_torus,
        base_point = list(

            auto = function(X, ...) {
                list(p = auto_base_point_torus(X), convergence = NA_integer_)
            },

            optimize = function(X, n_start = 10, maxit = 200, seed = 42, ...) {                
                optimize_base_point_torus(X=X,
                                          n_start=n_start,
                                          maxit=maxit,
                                          seed=seed)
            }
        )
    )
}

#' Check that data lies in the torus representation [0,1)^d
#'
#' Checks whether every value of \code{X} lies in \eqn{[0,1)}, the torus
#' representation expected by the torus functions of this package (see
#' \code{\link{manifold_torus}}).
#'
#' @param X A matrix or data frame of (nominally) torus data.
#'
#' @return Logical: \code{TRUE} if every value of \code{X} lies in
#'   \eqn{[0,1)}, \code{FALSE} otherwise.
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' is_in_torus(X)
#'
#' Y <- matrix(rnorm(20), nrow = 10, ncol = 2)  # not in [0,1)
#' is_in_torus(Y)
#'
#' @export
is_in_torus <- function(X) {

    return(all(X >= 0 & X < 1))

}

#' Map torus data to the tangent space at a base point
#'
#' Applies the logarithmic map to send torus data in \eqn{[0,1)^d}
#' into the tangent space \eqn{(-0.5, 0.5]^d} at base point \code{p}: for
#' each coordinate, \code{(x - p + 0.5) \%\% 1 - 0.5}, i.e. the shortest
#' signed offset from \code{p} on that coordinate's circle. Points on the
#' cut locus of \code{p} map to the boundary of \eqn{(-0.5, 0.5]^d}; see
#' \code{\link{is_valid_torus}}.
#'
#' @param X A matrix or data frame of torus data in \eqn{[0,1)^d}; rows are
#'   samples, columns are coordinates.
#' @param p Numeric vector of length \code{ncol(X)}: the base point on the
#'   torus.
#'
#' @return A numeric matrix of the same shape as \code{X}, with values in
#'   \eqn{(-0.5, 0.5]^d}.
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' log_map_torus(X, p = c(0.5, 0.5))
#'
#' @export
log_map_torus <- function(X, p) {

    X <- as.matrix(X)
    ((sweep(X, 2, p) + 0.5) %% 1) - 0.5

}

#' Torus distance from a data set to a point
#'
#' Computes the flat-torus distance between each row of \code{X} and the
#' point \code{p}: for each coordinate, the shortest wraparound distance
#' \code{min(|x-p|, 1-|x-p|)} on the unit circle, combined across
#' coordinates via the usual Euclidean norm.
#'
#' @param X A matrix or data frame of torus data in \eqn{[0,1)^d}; rows are
#'   samples, columns are coordinates.
#' @param p Numeric vector of length \code{ncol(X)}: a point on the torus.
#'
#' @return Numeric vector of length \code{nrow(X)}: the torus distance from
#'   each row of \code{X} to \code{p}.
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' dist_torus(X, c(0.5, 0.5))
#'
#' @export
dist_torus <- function(X, p) {

    X <- as.matrix(X)
    diff <- abs(sweep(X, 2, p))
    diff <- pmin(diff, 1 - diff)
    sqrt(rowSums(diff^2))
}

#' Check whether a base point is valid for a torus data set
#'
#' Checks that \code{p} has one coordinate per column of \code{X}, and
#' warns if any row of \code{X} lies on the cut locus of \code{p}, where
#' the inverse exponential map is undefined.
#'
#' @param X A matrix or data frame of torus data in \eqn{[0,1)^d}; rows are
#'   samples, columns are coordinates.
#' @param p Numeric vector of length \code{ncol(X)} giving the base point on
#'   the torus.
#'
#' @return Logical: \code{TRUE} if no row of \code{X} falls on the cut locus
#'   of \code{p}, \code{FALSE} otherwise. Stops with an error if \code{p}
#'   has the wrong length.
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' is_valid_torus(X, p = c(0, 0))
#'
#' @export
is_valid_torus <- function(X, p) {

    if (length(p) != ncol(X)) {
        stop(sprintf("The length of the base point (%d) mismatches the number of columns of X (%d)!", length(p), ncol(X)))
    }

    cut_locus_idx <- apply(X, 1, function(row) any(abs( (row - p) %% 1 - 0.5 )  < 1e-12))
    n_cut <- sum(cut_locus_idx)
    if (n_cut > 0) {
        warning(sprintf("%d samples detected on the cut locus. Consider modifying the base point.", n_cut))
    }

    invisible(n_cut == 0)
}

#' Midpoint of the largest gap on the circle [0,1)
#'
#' Given a set of forbidden values on the circle, returns the point that is
#' as far as possible from every forbidden value: the midpoint of the
#' largest circular gap between them. This guarantees the returned point is
#' not itself in the forbidden set.
#'
#' @param values Numeric vector in \eqn{[0,1)}.
#'
#' @return A scalar in \eqn{[0,1)} outside \code{values} (in the largest
#'   gap).
#'
#' @examples
#' largest_gap_point_torus(c(0.1, 0.3, 0.35))
#'
#' @export
largest_gap_point_torus <- function(values) {

    if (length(values) == 0) return(0.0)
    sv <- sort(values %% 1)
    n <- length(sv)

    gaps <- c(diff(sv), sv[1] + 1 - sv[n]) # Circular gaps between consecutive sorted values
    which_max <- which.max(gaps)

    if (which_max < n) {
        midpoint <- (sv[which_max] + sv[which_max + 1]) / 2
    } else {
        midpoint <- (sv[n] + sv[1] + 1) / 2
        midpoint <- midpoint %% 1
    }
    return(midpoint)
}

#' Find a torus base point via the maximal-gap method
#'
#' For the flat torus \eqn{T^d} with base point \code{p}, the cut locus is
#' \code{cut(p) = \{q : (q_j - p_j) mod 1 = 0.5 for some j\}}. A sample
#' point \code{q} is on the cut locus iff \code{p_j = (q_j - 0.5) mod 1}
#' for some coordinate \code{j}. To make the cut locus disjoint from the
#' sample, the function chooses, independently for each coordinate \code{j}, a
#' value \code{p_j} that is not equal to any of the forbidden values
#' \code{(q_ij - 0.5) mod 1}, i.e., a point in the largest circular gap
#' of the forbidden set. This yields a base point with zero cut-locus
#' overlap (up to floating-point tolerance) by construction, rather than
#' by search, but it is not centered on the data in any other sense.
#'
#' @param X A matrix or data frame of torus data in \eqn{[0,1)^d}; rows are
#'   samples, columns are coordinates.
#'
#' @return Numeric vector of length \code{ncol(X)}: the chosen base point.
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' auto_base_point_torus(X)
#'
#' @export
auto_base_point_torus <- function(X) {

    X <- as.matrix(X)
    apply(X, 2, function(col) largest_gap_point_torus((col - 0.5) %% 1))

}

#' Find a torus base point via multi-start optimization
#'
#' Torus case of \code{\link{optimize_base_point}}: numerically minimizes
#' the total torus distance (\code{\link{dist_torus}}) from a candidate
#' point \code{p} to every row of \code{X}, searching over the box
#' \eqn{[0,1]^d} with restarts drawn uniformly in it. See
#' \code{\link{choose_base_point}}.
#'
#' @param X A matrix or data frame of torus data in \eqn{[0,1)^d}; rows are
#'   samples, columns are coordinates.
#' @param n_start Number of random restarts used to avoid local minima of
#'   the (non-convex) objective. Defaults to 10.
#' @param maxit Maximum number of \code{\link[stats]{optim}} iterations per
#'   restart. Defaults to 200.
#' @param seed Integer random seed used to generate the restart points, for
#'   reproducibility. Defaults to 42.
#'
#' @return A list with three elements: \code{p} (the best point found
#'   across all restarts, clamped to \eqn{[0,1]^d}), \code{value} (its
#'   total torus distance to the rows of \code{X}), and \code{convergence}
#'   (the \code{\link[stats]{optim}} convergence code of the best restart;
#'   \code{0} indicates successful convergence).
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' optimize_base_point_torus(X, n_start = 5)
#'
#' @export
optimize_base_point_torus <- function(X, n_start = 10, maxit = 200, seed = 42) {

    d <- ncol(as.matrix(X))

    # Search over the closed box [0,1]^d rather than is_in_torus's [0,1)^d:
    # L-BFGS-B can land exactly on the upper bound 1 (the same point as 0)
    optimize_base_point(X = X,
                        dist = dist_torus,
                        init = function(n) matrix(stats::runif(n * d), ncol = d),
                        in_domain = function(p) all(p >= 0 & p <= 1),
                        lower = rep(0, d),
                        upper = rep(1, d),
                        project = function(p) pmin(pmax(p, 0), 1),
                        n_start = n_start,
                        maxit = maxit,
                        seed = seed)
}
