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
#' largest_gap_point(c(0.1, 0.3, 0.35))
#'
#' @export
largest_gap_point <- function(values) {

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

#' Find a base point via the maximal-gap method
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
#' auto_base_point(X)
#'
#' @export
auto_base_point <- function(X) {

    X <- as.matrix(X)
    apply(X, 2, function(col) largest_gap_point((col - 0.5) %% 1))

}

#' Find a base point via multi-start optimization
#'
#' Numerically minimizes \code{\link{torus_dist_sum}}, the total torus
#' distance from a candidate point \code{p} to every row of \code{X}, via
#' multi-start \code{\link[stats]{optim}} (L-BFGS-B, falling back to
#' Nelder-Mead if it errors). See \code{\link{choose_base_point}}.
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
#'   \code{\link{torus_dist_sum}} objective value), and \code{convergence}
#'   (the \code{\link[stats]{optim}} convergence code of the best restart;
#'   \code{0} indicates successful convergence).
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' optimize_base_point(X, n_start = 5)
#'
#' @export
optimize_base_point <- function(X, n_start = 10, maxit = 200, seed = 42) {

    X <- as.matrix(X)
    d <- ncol(X)

    # Box-constrained optimization using optim
    lower <- rep(0, d)
    upper <- rep(1, d)

    best_p <- NULL
    best_val <- Inf
    best_convergence <- NA_integer_

    # Generate initial points: uniformly distributed in [0,1]^d
    set.seed(seed)
    init_points <- matrix(stats::runif(n_start * d), ncol = d)

    for (k in seq_len(n_start)) {

        start <- init_points[k, ]
        # Use the L-BFGS-B method to handle boundaries; fall back to Nelder-Mead if unavailable
        opt <- tryCatch(
            stats::optim(par = start,
                fn = torus_dist_sum,
                X = X,
                method = "L-BFGS-B",
                lower = lower,
                upper = upper,
                control = list(maxit = maxit)),
            error = function(e) {
                stats::optim(par = start,
                    fn = torus_dist_sum,
                    X = X,
                    method = "Nelder-Mead",
                    control = list(maxit = maxit))
            }
        )
        if (opt$value < best_val) {
            best_val <- opt$value
            best_p <- opt$par
            best_convergence <- opt$convergence
        }
    }

    # Ensure the optimal solution lies within [0,1]
    best_p <- pmin(pmax(best_p, 0), 1)

    list(p = best_p, value = best_val, convergence = best_convergence)
}

#' Choose a base point for the logarithmic map
#'
#' Finds a base point \code{p} on the torus \eqn{[0,1)^d} for \code{X},
#' using one of two methods:
#' \itemize{
#'   \item \code{"auto"} (the default): the closed-form, per-coordinate
#'     maximal-gap method (\code{\link{auto_base_point}}), which is
#'     \emph{not} centered on the data, but is guaranteed to avoid the cut
#'     locus of \code{X} by construction, and is much cheaper to compute.
#'   \item \code{"optimize"}: numerically minimizes
#'     \code{\link{torus_dist_sum}}, the total torus distance from
#'     \code{p} to every row of \code{X} (a geometric-median-style,
#'     data-centered point), via multi-start \code{\link[stats]{optim}}.
#'     Cut-locus avoidance is checked afterward, not guaranteed by
#'     construction.
#' }
#' Either way, the chosen point is validated with \code{\link{is_valid}}
#' before being returned.
#'
#' @param X A matrix or data frame of torus data in \eqn{[0,1)^d}; rows are
#'   samples, columns are coordinates.
#' @param method Either \code{"auto"} (the default) or \code{"optimize"};
#'   see above.
#' @param n_start Number of random restarts for \code{method = "optimize"},
#'   used to avoid local minima of the (non-convex) objective. Defaults to
#'   10. Ignored for \code{method = "auto"}.
#' @param maxit Maximum number of \code{\link[stats]{optim}} iterations per
#'   restart, for \code{method = "optimize"}. Defaults to 200. Ignored for
#'   \code{method = "auto"}.
#' @param seed Integer random seed used to generate the restart points, for
#'   \code{method = "optimize"}. Defaults to 42. Ignored for
#'   \code{method = "auto"}.
#'
#' @return A list with three elements:
#'   \item{p}{Numeric vector of length \code{ncol(X)}: the chosen base
#'     point, clamped to \eqn{[0,1]^d}.}
#'   \item{value}{The total torus distance from \code{p} to every row of
#'     \code{X} (\code{\link{torus_dist_sum}}), computed for either method
#'     so the two are comparable. Lower means \code{p} is more central.}
#'   \item{convergence}{For \code{method = "optimize"}, the
#'     \code{\link[stats]{optim}} convergence code of the best restart
#'     (\code{0} indicates successful convergence). \code{NA} for
#'     \code{method = "auto"}, which involves no iterative optimization.}
#'
#'   Stops with an error (checked via \code{\link{is_valid}}) if the chosen
#'   \code{p} happens to land on the cut locus of some row of \code{X}: for
#'   \code{method = "optimize"} this suggests a higher \code{n_start},
#'   while for \code{method = "auto"} it should not happen (up to
#'   floating-point tolerance) and indicates degenerate data.
#'
#' @examples
#' X <- matrix(stats::runif(20), nrow = 10, ncol = 2)
#'
#' result_auto <- choose_base_point(X)
#' result_auto$p
#'
#' result <- choose_base_point(X, method = "optimize", n_start = 5)
#' result$p
#'
#' @export
choose_base_point <- function(X, method = c("auto", "optimize"),
                               n_start = 10, maxit = 200, seed = 42) {

    method <- match.arg(method)
    X <- as.matrix(X)
    d <- ncol(X)
    if (d == 0) stop("X must have at least 1 column.")

    if (method == "auto") {

        best_p <- auto_base_point(X)
        best_val <- torus_dist_sum(best_p, X)
        best_convergence <- NA_integer_

    } else {

        opt_result <- optimize_base_point(X=X,
                                        n_start=n_start,
                                        maxit=maxit,
                                        seed=seed)
        best_p <- opt_result$p
        best_val <- opt_result$value
        best_convergence <- opt_result$convergence
    }

    best_valid <- is_valid(X, best_p)
    if (!best_valid) {
        if (method == "optimize") {
            stop("The chosen base point lies on the cut locus of some sample point. Consider increasing n_start and retrying.")
        } else {
            stop("The chosen base point lies on the cut locus of some sample point despite the maximal-gap guarantee; the data may be degenerate.")
        }
    }

    return(list(p = best_p, value = best_val, convergence = best_convergence))
}
