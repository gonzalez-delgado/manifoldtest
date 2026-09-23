#' Distance covariance independence test for torus-valued data.
#'
#' Aligns \code{X} and \code{Y} to a common sample size, chooses (or
#' validates a user-supplied) base point for each, maps both data sets from
#' the torus to the tangent space via the inverse exponential map, then runs
#' a distance covariance permutation test to assess independence between
#' them.
#'
#' @param X A matrix or data frame of torus data in \eqn{[0,1)^d}; rows
#'   are samples, columns are coordinates.
#' @param Y A matrix or data frame of torus data in \eqn{[0,1)^d}; rows
#'   are samples, columns are coordinates.
#' @param px_base Numeric vector of length \code{ncol(X)} giving the base
#'   point used to map \code{X} to the tangent space. If \code{NULL} (the
#'   default) or not \code{\link{is_valid}} for \code{X} (i.e. some row of
#'   \code{X} lies on its cut locus), a base point is chosen automatically
#'   with \code{\link{choose_base_point}}.
#' @param py_base Numeric vector of length \code{ncol(Y)} giving the base
#'   point used to map \code{Y} to the tangent space. If \code{NULL} (the
#'   default) or not \code{\link{is_valid}} for \code{Y}, a base point is
#'   chosen automatically with \code{\link{choose_base_point}}.
#' @param R Number of permutations for the distance covariance test, passed
#'   to \code{\link[energy]{dcov.test}}. Defaults to 500.
#' @param seed Random seed used by \code{\link{align_sample_sizes}} when the
#'   sample sizes of \code{X} and \code{Y} differ, by
#'   \code{\link{choose_base_point}} when \code{method = "optimize"}, and by
#'   \code{\link{center_outward_transform}} when \code{center_outward =
#'   TRUE}. Defaults to \code{2024}.
#' @param center_outward If \code{TRUE} (the default), applies the
#'   center-outward rank transform (Hallin et al. 2021,
#'   \code{\link{center_outward_transform}}) to the log-mapped \code{X} and
#'   \code{Y} before testing, so the distance covariance test is computed
#'   on multidimensional ranks bounded in the unit ball rather than on the
#'   raw tangent-space coordinates. This makes the test valid without
#'   moment conditions on the original data. If \code{FALSE}, the test runs
#'   directly on the tangent-space coordinates from \code{\link{log_map}}.
#' @param method Method passed to \code{\link{choose_base_point}} for any
#'   base point that must be chosen automatically: \code{"auto"} (the
#'   default, closed-form maximal-gap method) or \code{"optimize"}
#'   (numerical, data-centered). See \code{\link{choose_base_point}}.
#' @param n_start Number of random restarts passed to
#'   \code{\link{choose_base_point}}. Only used when \code{method =
#'   "optimize"}. Defaults to 10.
#' @param maxit Maximum \code{\link[stats]{optim}} iterations per restart,
#'   passed to \code{\link{choose_base_point}}. Only used when
#'   \code{method = "optimize"}. Defaults to 200.
#' @param verbose If \code{TRUE} (the default), prints a summary of the
#'   alignment, base points, and test result to the console.
#'
#' @return A list with five elements:
#'   \item{dcov_test}{The \code{dcov.test} result object, as returned by
#'     \code{\link[energy]{dcov.test}}.}
#'   \item{px_base}{The base point used to map \code{X}: either the
#'     supplied \code{px_base}, or the one chosen automatically.}
#'   \item{py_base}{The base point used to map \code{Y}, analogous to
#'     \code{px_base}.}
#'   \item{X_mapped}{The mapped \code{X} samples actually used for the
#'     test: rows in the unit ball \eqn{S^d} from
#'     \code{\link{center_outward_transform}} if \code{center_outward =
#'     TRUE}, otherwise the tangent-space coordinates from
#'     \code{\link{log_map}}.}
#'   \item{Y_mapped}{The mapped \code{Y} samples actually used for the
#'     test, analogous to \code{X_mapped}.}
#'
#' @examples
#' \donttest{
#' X <- matrix(runif(30), nrow = 10, ncol = 3)
#' Y <- matrix(runif(20), nrow = 10, ncol = 2)
#' result <- torus.dcov.test(X, Y, R = 200)
#' result$dcov_test$p.value
#' result$px_base
#' dim(result$X_mapped)
#' }
#'
#' @export
torus.dcov.test <- function(X, Y, px_base = NULL, py_base = NULL,
                            R = 500, seed = 2024, center_outward = TRUE,
                            method = c("auto", "optimize"), n_start = 10, maxit = 200,
                            verbose = TRUE) {

    method <- match.arg(method)

    # -------------- Input validation --------------

    if (!is_in_torus(X)) stop("X should lie in [0,1)^d (torus representation).")
    if (!is_in_torus(Y)) stop("Y should lie in [0,1)^d (torus representation).")

    # -------------- Sample size alignment --------------

    aligned <- align_sample_sizes(X_raw=X,
                                Y_raw=Y,
                                seed=seed)
    X <- aligned$X
    Y <- aligned$Y

    # -------------- Base point selection --------------

    if (is.null(px_base) || !is_valid(X, px_base)) {

        if(verbose){message("Choosing base point for X...\n")}
        px_base <- choose_base_point(X=X,
                                    method=method,
                                    n_start=n_start,
                                    maxit=maxit,
                                    seed=seed)$p
        if(verbose){message(sprintf("Chosen base point for X: (%s)\n", paste(sprintf("%.2f", px_base), collapse = ", ")))}
    }
    if (is.null(py_base) || !is_valid(Y, py_base)) {

        if(verbose){message("Choosing base point for Y...\n")}
        py_base <- choose_base_point(X=Y,
                                    method=method,
                                    n_start=n_start,
                                    maxit=maxit,
                                    seed=seed)$p
        if(verbose){message(sprintf("Chosen base point for Y: (%s)\n", paste(sprintf("%.2f", py_base), collapse = ", ")))}
    }

    # -------------- Mapping to tangent space --------------

    X_euclidean <- log_map(X, px_base)
    Y_euclidean <- log_map(Y, py_base)

    if(center_outward){

        # -------------- Center-outward transform --------------

        if(verbose){message("Performing center-outward transform...\n")}
        X_test <- center_outward_transform(X_euclidean, seed=seed)
        Y_test <- center_outward_transform(Y_euclidean, seed=seed)

    }else{

        X_test <- X_euclidean
        Y_test <- Y_euclidean
    }

    # -------------- Distance Covariance Test --------------

    if(verbose){message("Performing distance covariance test...\n")}
    test_result <- energy::dcov.test(X_test, Y_test, R = R)

    if (verbose) {
        message("\n========================================")
        message(" Distance Covariance Independence Test Result")
        message("========================================")
        message(sprintf("Sample size: %d\n", nrow(X_euclidean)))
        message(sprintf("X Dimension: %d\n", ncol(X_euclidean)))
        message(sprintf("Sample size: %d\n", nrow(Y_euclidean)))
        message(sprintf("Y Dimension: %d\n", ncol(Y_euclidean)))
        message(sprintf("Distance covariance statistic V_n: %.5f", test_result$statistic))
        message(sprintf("p: %.5f\n", test_result$p.value))
        message("========================================")
    }

    return(list(
        dcov_test = test_result,
        px_base   = px_base,
        py_base   = py_base,
        X_mapped  = X_test,
        Y_mapped  = Y_test
    ))
}
