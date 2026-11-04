

#' Align sample sizes of two data sets
#'
#' Checks whether \code{X_raw} and \code{Y_raw} have the same number of rows,
#' and if not, randomly subsamples both to the smaller of the two sizes so
#' that they can be compared or paired downstream.
#'
#' @param X_raw A matrix or data frame; rows are samples.
#' @param Y_raw A matrix or data frame; rows are samples.
#' @param seed Integer random seed used when subsampling, for reproducibility.
#'   Defaults to \code{2024}.
#'
#' @return A list with two elements, \code{X} and \code{Y}, each a matrix or
#'   data frame with matching row counts. If the inputs already had equal
#'   row counts, they are returned unchanged.
#'
#' @examples
#' X_raw <- matrix(rnorm(30), nrow = 10)
#' Y_raw <- matrix(rnorm(24), nrow = 8)
#' aligned <- align_sample_sizes(X_raw, Y_raw)
#' nrow(aligned$X)
#' nrow(aligned$Y)
#'
#' @export
align_sample_sizes <- function(X_raw, Y_raw, seed = 2024) {

    n_X <- nrow(X_raw)
    n_Y <- nrow(Y_raw)
    if (n_X != n_Y) {
    
        message(sprintf("Warning: The sample size of X (%d) is inconsistent with that of Y (%d), and will be automatically aligned to the smaller value. \n", n_X, n_Y))
        n_min <- min(n_X, n_Y)
        set.seed(seed)
    
        # Randomly remove redundant samples
        idx_X <- sample(1:n_X, n_min)
        idx_Y <- sample(1:n_Y, n_min)
        X <- X_raw[idx_X, , drop = FALSE]
        Y <- Y_raw[idx_Y, , drop = FALSE]
        message(sprintf("Sample size after alignment: %d\n", n_min))
    } else {
        X <- X_raw
        Y <- Y_raw
    }
    return(list(X = X, Y = Y))
}

#' Find a base point via multi-start optimization
#'
#' Manifold-agnostic search for a data-centered base point: numerically
#' minimizes the total geodesic distance \code{sum(dist(X, p))} from a
#' candidate point \code{p} to every row of \code{X}, via multi-start
#' \code{\link[stats]{optim}} (L-BFGS-B, falling back to Nelder-Mead if it
#' errors). The manifold enters only through \code{dist} and through the
#' search space, which is described by \code{init}, \code{in_domain},
#' \code{lower}, \code{upper} and \code{project}. See
#' \code{\link{optimize_base_point_torus}} for the torus case.
#'
#' @param X A matrix or data frame of manifold data; rows are samples,
#'   columns are coordinates.
#' @param dist Function \code{f(X, p)} returning the vector of geodesic
#'   distances between each row of \code{X} and \code{p} (see
#'   \code{\link{new_manifold}}).
#' @param init Function \code{f(n)} returning a matrix whose \code{n} rows
#'   are the restart points. It is called after \code{set.seed(seed)}.
#' @param in_domain Function \code{f(p)} returning \code{TRUE} if \code{p}
#'   is an admissible candidate; the objective is \code{Inf} otherwise.
#'   Defaults to accepting every \code{p}.
#' @param lower,upper Bounds on \code{p} passed to
#'   \code{\link[stats]{optim}} for L-BFGS-B. Default to \code{-Inf} and
#'   \code{Inf} (unconstrained).
#' @param project Function \code{f(p)} applied to the best point found,
#'   e.g. to map it back onto the manifold. Defaults to
#'   \code{\link[base]{identity}}.
#' @param n_start Number of random restarts used to avoid local minima of
#'   the (non-convex) objective. Defaults to 10.
#' @param maxit Maximum number of \code{\link[stats]{optim}} iterations per
#'   restart. Defaults to 200.
#' @param seed Integer random seed used to generate the restart points, for
#'   reproducibility. Defaults to 42.
#'
#' @return A list with three elements: \code{p} (the best point found
#'   across all restarts, after \code{project}), \code{value} (its total
#'   distance to the rows of \code{X}), and \code{convergence} (the
#'   \code{\link[stats]{optim}} convergence code of the best restart;
#'   \code{0} indicates successful convergence).
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' optimize_base_point(X, dist = dist_torus,
#'                     init = function(n) matrix(runif(2 * n), ncol = 2),
#'                     lower = c(0, 0), upper = c(1, 1), n_start = 5)
#'
#' @export
optimize_base_point <- function(X, dist, init, in_domain = function(p) TRUE,
                                lower = -Inf, upper = Inf, project = identity,
                                n_start = 10, maxit = 200, seed = 42) {

    X <- as.matrix(X)

    # Objective: total distance to the data, infinite outside the domain
    objective <- function(p) {
        if (!in_domain(p)) return(Inf)
        sum(dist(X, p))
    }

    best_p <- NULL
    best_val <- Inf
    best_convergence <- NA_integer_

    set.seed(seed)
    init_points <- as.matrix(init(n_start))
    if (nrow(init_points) != n_start) {
        stop(sprintf("init(%d) returned %d restart points instead of %d.", n_start, nrow(init_points), n_start))
    }

    for (k in seq_len(n_start)) {

        start <- init_points[k, ]
        # Use the L-BFGS-B method to handle bounds; fall back to Nelder-Mead if it errors
        opt <- tryCatch(
            stats::optim(par = start,
                fn = objective,
                method = "L-BFGS-B",
                lower = lower,
                upper = upper,
                control = list(maxit = maxit)),
            error = function(e) {
                stats::optim(par = start,
                    fn = objective,
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

    list(p = project(best_p), value = best_val, convergence = best_convergence)
}
