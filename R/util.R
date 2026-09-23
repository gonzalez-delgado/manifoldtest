

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


#' Check whether a base point is valid for a torus data set
#'
#' Checks that \code{p_base} has one coordinate per column of \code{X}, and
#' warns if any row of \code{X} lies on the cut locus of \code{p_base}, where
#' the inverse exponential map is undefined.
#'
#' @param X A matrix or data frame of torus data in \eqn{[0,1)^d}; rows are
#'   samples, columns are coordinates.
#' @param p_base Numeric vector of length \code{ncol(X)} giving the base
#'   point on the torus.
#'
#' @return Logical: \code{TRUE} if no row of \code{X} falls on the cut locus
#'   of \code{p_base}, \code{FALSE} otherwise. Stops with an error if
#'   \code{p_base} has the wrong length. As a side effect, emits a warning if
#'   the point is invalid, or a message if it is valid.
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' is_valid(X, p_base = c(0, 0))
#'
#' @export
is_valid <- function(X, p_base) {

    if (length(p_base) != ncol(X)) {
        stop(sprintf("The length of p_base (%d) mismatches the number of columns of X (%d)!", length(p_base), ncol(X)))
    }

    cut_locus_idx <- apply(X, 1, function(row) any(abs( (row - p_base) %% 1 - 0.5 )  < 1e-12))
    n_cut <- sum(cut_locus_idx)
    if (n_cut > 0) {
        warning(sprintf("%d samples detected on the cut locus. Consider modifying the base point.", n_cut))
    } else {
        message("No sample detected on the cut locus.")
    }

    invisible(n_cut == 0)
}

#' Total torus distance from a candidate point to a data set
#'
#' Computes the sum, over every row of \code{X}, of the flat-torus distance
#' between that row and the candidate point \code{p}: for each coordinate,
#' the shortest wraparound distance \code{min(|x-p|, 1-|x-p|)} on the unit
#' circle, combined across coordinates via the usual Euclidean norm.
#'
#' @param p Numeric vector of length \code{ncol(X)}: the candidate point on
#'   the torus \eqn{[0,1)^d}.
#' @param X A numeric matrix of torus data in \eqn{[0,1)^d}; rows are
#'   samples, columns are coordinates.
#'
#' @return A single number: the total torus distance from \code{p} to every
#'   row of \code{X}, or \code{Inf} if any coordinate of \code{p} falls
#'   outside \eqn{[0,1]}.
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' torus_dist_sum(c(0.5, 0.5), X)
#'
#' @export
torus_dist_sum <- function(p, X) {

    if (any(p < 0 | p > 1)) return(Inf)
    diff <- abs(sweep(X, 2, p))
    diff <- pmin(diff, 1 - diff)
    sum(sqrt(rowSums(diff^2)))
}

#' Check that data lies in the torus representation [0,1)^d
#'
#' Checks whether every value of \code{X} lies in \eqn{[0,1)}, the torus
#' representation expected throughout this package (see
#' \code{\link{log_map}}, \code{\link{is_valid}},
#' \code{\link{torus_dist_sum}}).
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
#' \code{\link{is_valid}}.
#'
#' @param X_torus A matrix or data frame of torus data in \eqn{[0,1)^d};
#'   rows are samples, columns are coordinates.
#' @param p Numeric vector of length \code{ncol(X_torus)}: the base point on
#'   the torus.
#'
#' @return A numeric matrix of the same shape as \code{X_torus}, with values
#'   in \eqn{(-0.5, 0.5]^d}.
#'
#' @examples
#' X <- matrix(runif(20), nrow = 10, ncol = 2)
#' log_map(X, p = c(0.5, 0.5))
#'
#' @export
log_map <- function(X_torus, p) {

    X_torus <- as.matrix(X_torus)
    ((sweep(X_torus, 2, p) + 0.5) %% 1) - 0.5

}

