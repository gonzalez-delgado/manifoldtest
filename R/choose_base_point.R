#' Choose a base point for the logarithmic map
#'
#' Finds a base point \code{p} on the manifold for \code{X}, using one of
#' the base point methods provided by the manifold (see
#' \code{\link{new_manifold}}), then validates it with the manifold's
#' cut-locus check (\code{is_valid}). For the torus
#' (\code{\link{manifold_torus}}) the methods are:
#' \itemize{
#'   \item \code{"auto"} (the default): the closed-form, per-coordinate
#'     maximal-gap method (\code{\link{auto_base_point_torus}}), which is
#'     \emph{not} centered on the data, but is guaranteed to avoid the cut
#'     locus of \code{X} by construction, and is much cheaper to compute.
#'   \item \code{"optimize"}: numerically minimizes the total torus
#'     distance from \code{p} to every row of \code{X} (a
#'     geometric-median-style, data-centered point), via multi-start
#'     \code{\link[stats]{optim}} (\code{\link{optimize_base_point_torus}}).
#'     Cut-locus avoidance is checked afterward, not guaranteed by
#'     construction.
#' }
#'
#' @param X A matrix or data frame of manifold data; rows are samples,
#'   columns are coordinates.
#' @param manifold The manifold \code{X} lies on: a name accepted by
#'   \code{\link{get_manifold}} (currently only \code{"torus"}, the
#'   default) or an object created with \code{\link{new_manifold}}.
#' @param method Name of one of the manifold's base point methods. If
#'   \code{NULL} (the default), the manifold's first (default) method is
#'   used, i.e. \code{"auto"} for the torus.
#' @param n_start Number of random restarts for iterative methods (e.g.
#'   \code{"optimize"}), used to avoid local minima of the objective.
#'   Defaults to 10. Ignored by methods that do not use it.
#' @param maxit Maximum number of \code{\link[stats]{optim}} iterations per
#'   restart, for iterative methods. Defaults to 200. Ignored by methods
#'   that do not use it.
#' @param seed Integer random seed for iterative methods. Defaults to 42.
#'   Ignored by methods that do not use it.
#'
#' @return A list with three elements:
#'   \item{p}{The chosen base point.}
#'   \item{value}{The total geodesic distance from \code{p} to every row of
#'     \code{X}, computed with the manifold's \code{dist} for every method
#'     so they are comparable. Lower means \code{p} is more central.}
#'   \item{convergence}{For iterative methods, the
#'     \code{\link[stats]{optim}} convergence code of the best restart
#'     (\code{0} indicates successful convergence). \code{NA} for
#'     non-iterative methods such as \code{"auto"}.}
#'
#'   Stops with an error if the chosen \code{p} happens to land on the cut
#'   locus of some row of \code{X}: for \code{method = "optimize"} this
#'   suggests a higher \code{n_start}, while for the torus \code{method =
#'   "auto"} it should not happen (up to floating-point tolerance) and
#'   indicates degenerate data.
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
choose_base_point <- function(X, manifold = "torus", method = NULL,
                               n_start = 10, maxit = 200, seed = 42) {

    manifold <- get_manifold(manifold)
    if (is.null(method)) method <- names(manifold$base_point)[1]
    method <- match.arg(method, names(manifold$base_point))

    X <- as.matrix(X)
    d <- ncol(X)
    if (d == 0) stop("X must have at least 1 column.")

    result <- manifold$base_point[[method]](X = X,
                                            n_start = n_start,
                                            maxit = maxit,
                                            seed = seed)
    best_p <- result$p
    best_val <- sum(manifold$dist(X, best_p))
    best_convergence <- result$convergence

    best_valid <- manifold$is_valid(X, best_p)
    if (!best_valid) {
        if (method == "optimize") {
            stop("The chosen base point lies on the cut locus of some sample point. Consider increasing n_start and retrying.")
        } else {
            stop(sprintf("The base point chosen with method '%s' lies on the cut locus of some sample point; the data may be degenerate.", method))
        }
    }

    return(list(p = best_p, value = best_val, convergence = best_convergence))
}
