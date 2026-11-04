#' The Grassmann manifold.
#'
#' Returns the \code{\link{new_manifold}} object describing the Grassmann
#' manifold \eqn{Gr(n, k)} of \code{k}-dimensional linear subspaces of
#' \eqn{R^n}. A subspace is represented by an \code{n} x \code{k} matrix
#' \eqn{U} with orthonormal columns spanning it (a point of the Stiefel
#' manifold, \eqn{U^T U = I_k}), stored as a row of length \code{n * k}
#' holding \code{as.vector(U)}, i.e. \eqn{U} column by column. Data on
#' \eqn{Gr(n, k)} are thus matrices with \code{n * k} columns, and base
#' points are vectors of length \code{n * k}. Two matrices \eqn{U} and
#' \eqn{UR}, with \eqn{R} orthogonal, represent the same subspace; all
#' functions are invariant to this choice up to an isometry of the tangent
#' space, which leaves distance covariance unchanged.
#'
#' Since \code{n} and \code{k} cannot be inferred from the data, the
#' manifold must be passed as an object, e.g. \code{manifold =
#' manifold_grassmann(4, 2)}, rather than by name.
#'
#' The tangent space at a base point is identified with
#' \eqn{R^{(n-k)k}} (see \code{\link{log_map_grassmann}}) and the cut locus
#' of \eqn{U_0} is the set of subspaces with at least one principal angle
#' equal to \eqn{\pi/2} with it. Its components are
#' \code{\link{is_in_grassmann}}, \code{\link{log_map_grassmann}},
#' \code{\link{dist_grassmann}}, \code{\link{is_valid_grassmann}}, and two
#' base point methods: \code{"auto"} (the default, the closed-form
#' projection mean \code{\link{auto_base_point_grassmann}}) and
#' \code{"optimize"} (\code{\link{optimize_base_point_grassmann}}).
#'
#' @param n Integer: dimension of the ambient space.
#' @param k Integer: dimension of the subspaces, \code{1 <= k <= n - 1}.
#'
#' @return An object of class \code{"manifold"}.
#'
#' @examples
#' manifold_grassmann(4, 2)
#'
#' @export
manifold_grassmann <- function(n, k) {

    check_grassmann_dims(n, k)
    force(n); force(k)

    new_manifold(
        name = sprintf("Grassmann manifold Gr(%d, %d)", n, k),
        representation = sprintf("rows as.vector(U) of %d x %d matrices U with orthonormal columns", n, k),
        contains = function(X) is_in_grassmann(X, n, k),
        log_map = function(X, p) log_map_grassmann(X, p, n, k),
        dist = function(X, p) dist_grassmann(X, p, n, k),
        is_valid = function(X, p) is_valid_grassmann(X, p, n, k),
        base_point = list(

            auto = function(X, ...) {
                list(p = auto_base_point_grassmann(X, n, k), convergence = NA_integer_)
            },

            optimize = function(X, n_start = 10, maxit = 200, seed = 42, ...) {
                optimize_base_point_grassmann(X=X,
                                              n=n,
                                              k=k,
                                              n_start=n_start,
                                              maxit=maxit,
                                              seed=seed)
            }
        )
    )
}

check_grassmann_dims <- function(n, k) {

    if (!is.numeric(n) || !is.numeric(k) || length(n) != 1 || length(k) != 1 ||
        n != round(n) || k != round(k) || k < 1 || k > n - 1) {
        stop("n and k must be integers with 1 <= k <= n - 1.")
    }
}

# Reshape a base point given as a vector of length n * k (or an n x k
# matrix) into an n x k matrix, checking its shape and orthonormality.
as_stiefel <- function(p, n, k, tol = 1e-6) {

    if (length(p) != n * k) {
        stop(sprintf("The length of the base point (%d) mismatches n * k = %d!", length(p), n * k))
    }
    U <- matrix(p, n, k)
    if (max(abs(crossprod(U) - diag(k))) > tol) {
        stop("The base point must have orthonormal columns (t(U) %*% U = I).")
    }
    U
}

# Orthonormal basis of the orthogonal complement of span(U0), used to
# identify the tangent space at U0 with R^((n-k)k). Deterministic in U0.
grassmann_complement <- function(U0) {

    k <- ncol(U0)
    qr.Q(qr(U0), complete = TRUE)[, -seq_len(k), drop = FALSE]
}

# Riemannian logarithm at U0 of one Stiefel matrix V, as the n x k
# horizontal lift Delta (Bendokat, Zimmermann and Absil, 2021, Alg. 5.3).
grassmann_log_one <- function(V, U0) {

    k <- ncol(U0)
    s1 <- svd(crossprod(V, U0))
    Vs <- V %*% s1$u %*% t(s1$v)                 # Procrustes-aligned representative
    s2 <- svd(Vs - U0 %*% crossprod(U0, Vs))
    s2$u %*% diag(asin(pmin(s2$d, 1)), k, k) %*% t(s2$v)
}

#' Check that data lies on the Grassmann manifold.
#'
#' Checks that \code{X} has \code{n * k} columns and that every row,
#' reshaped as an \code{n} x \code{k} matrix \eqn{U} (column by column),
#' has orthonormal columns, \eqn{U^T U = I_k}, up to \code{tol} (see
#' \code{\link{manifold_grassmann}}).
#'
#' @param X A matrix or data frame of (nominally) Grassmann data; rows are
#'   samples, with \code{n * k} columns.
#' @inheritParams manifold_grassmann
#' @param tol Tolerance on the entries of \eqn{U^T U - I_k}. Defaults to
#'   \code{1e-6}.
#'
#' @return Logical: \code{TRUE} if every row of \code{X} represents a point
#'   of \eqn{Gr(n, k)}, \code{FALSE} otherwise.
#'
#' @examples
#' X <- exp_map_grassmann(matrix(rnorm(40, sd = 0.3), nrow = 10), diag(4)[, 1:2])
#' is_in_grassmann(X, n = 4, k = 2)
#'
#' is_in_grassmann(matrix(runif(80), nrow = 10), n = 4, k = 2)
#'
#' @export
is_in_grassmann <- function(X, n, k, tol = 1e-6) {

    check_grassmann_dims(n, k)
    X <- as.matrix(X)
    if (ncol(X) != n * k || !is.numeric(X)) return(FALSE)
    all(apply(X, 1, function(row) {
        U <- matrix(row, n, k)
        max(abs(crossprod(U) - diag(k))) <= tol
    }))
}

#' Map Grassmann data to the tangent space at a base point.
#'
#' Applies the Riemannian logarithm of \eqn{Gr(n, k)} at base point
#' \eqn{U_0}, computed with Algorithm 5.3 of Bendokat, Zimmermann and
#' Absil (2021). The result for each row is the horizontal lift
#' \eqn{\Delta}, an \code{n} x \code{k} matrix with \eqn{U_0^T \Delta = 0}
#' whose singular values are the principal angles between the two
#' subspaces. It is expressed in the coordinates
#' \eqn{B = U_\perp^T \Delta}, where \eqn{U_\perp} is a fixed orthonormal
#' basis of the orthogonal complement of \eqn{U_0}, and returned as
#' \code{as.vector(B)}, of length \eqn{(n-k)k}. This identification is an
#' isometry: the Euclidean norm of the coordinates is the geodesic distance
#' to \eqn{U_0}. Points on the cut locus of \eqn{U_0} (see
#' \code{\link{is_valid_grassmann}}) are mapped to the boundary of the cut
#' domain, where the logarithm is not unique.
#'
#' @param X A matrix or data frame of Grassmann data; rows are samples,
#'   with \code{n * k} columns (see \code{\link{manifold_grassmann}}).
#' @param p The base point: a vector of length \code{n * k} or an \code{n}
#'   x \code{k} matrix with orthonormal columns.
#' @inheritParams manifold_grassmann
#'
#' @return A numeric matrix with \code{nrow(X)} rows and \code{(n - k) * k}
#'   columns: the tangent-space coordinates of each row of \code{X}.
#'
#' @references Bendokat, T., Zimmermann, R. and Absil, P.-A. (2021). A
#'   Grassmann manifold handbook: basic geometry and computational aspects.
#'   arXiv:2011.13699.
#'
#' @examples
#' U0 <- diag(4)[, 1:2]
#' X <- exp_map_grassmann(matrix(rnorm(40, sd = 0.3), nrow = 10), U0)
#' log_map_grassmann(X, U0, n = 4, k = 2)
#'
#' @export
log_map_grassmann <- function(X, p, n, k) {

    check_grassmann_dims(n, k)
    U0 <- as_stiefel(p, n, k)
    Up <- grassmann_complement(U0)
    X <- as.matrix(X)
    out <- t(apply(X, 1, function(row) {
        as.vector(crossprod(Up, grassmann_log_one(matrix(row, n, k), U0)))
    }))
    # apply() drops to a vector when (n - k) * k == 1
    matrix(out, nrow = nrow(X), ncol = (n - k) * k)
}

#' Map tangent vectors to the Grassmann manifold.
#'
#' Exponential map of \eqn{Gr(n, k)} at base point \eqn{U_0}, the inverse
#' of \code{\link{log_map_grassmann}} on the cut domain. Each row of
#' \code{V} holds tangent coordinates \code{as.vector(B)}, with \code{B} of
#' size \code{(n - k)} x \code{k}; the horizontal lift is
#' \eqn{\Delta = U_\perp B} and, writing its compact singular value
#' decomposition \eqn{\Delta = Q \Sigma R^T}, the point is
#' \eqn{U_0 R \cos(\Sigma) R^T + Q \sin(\Sigma) R^T}. Useful to simulate
#' data on \eqn{Gr(n, k)}.
#'
#' @param V Numeric matrix with \code{(n - k) * k} columns: tangent-space
#'   coordinates at \code{p}, one row per point.
#' @param p The base point: an \code{n} x \code{k} matrix with orthonormal
#'   columns, or a vector of length \code{n * k} (then \code{n} and
#'   \code{k} must be given).
#' @param n,k Integers: dimensions of \eqn{Gr(n, k)}. Default to the
#'   dimensions of \code{p} when it is a matrix.
#'
#' @return A numeric matrix with \code{nrow(V)} rows and \code{n * k}
#'   columns: Stiefel representatives \code{as.vector(U)} of the points.
#'
#' @examples
#' X <- exp_map_grassmann(matrix(rnorm(40, sd = 0.3), nrow = 10), diag(4)[, 1:2])
#' is_in_grassmann(X, n = 4, k = 2)
#'
#' @export
exp_map_grassmann <- function(V, p, n = nrow(p), k = ncol(p)) {

    if (is.null(n) || is.null(k)) stop("n and k must be given when p is a vector.")
    check_grassmann_dims(n, k)
    U0 <- as_stiefel(p, n, k)
    Up <- grassmann_complement(U0)
    V <- as.matrix(V)
    if (ncol(V) != (n - k) * k) {
        stop(sprintf("V must have (n - k) * k = %d columns.", (n - k) * k))
    }
    out <- t(apply(V, 1, function(row) {
        s <- svd(Up %*% matrix(row, n - k, k))
        U <- U0 %*% s$v %*% diag(cos(s$d), k, k) %*% t(s$v) +
            s$u %*% diag(sin(s$d), k, k) %*% t(s$v)
        as.vector(U)
    }))
    matrix(out, nrow = nrow(V), ncol = n * k)
}

#' Grassmann distance from a data set to a point.
#'
#' Computes the geodesic distance on \eqn{Gr(n, k)} between each row of
#' \code{X} and the base point \code{p}: the square root of the sum of the
#' squared principal angles between the two subspaces. It is evaluated as
#' the norm of \code{\link{log_map_grassmann}}, which is numerically
#' accurate for nearby subspaces.
#'
#' @inheritParams log_map_grassmann
#'
#' @return Numeric vector of length \code{nrow(X)}: the Grassmann distance
#'   from each row of \code{X} to \code{p}.
#'
#' @examples
#' U0 <- diag(4)[, 1:2]
#' X <- exp_map_grassmann(matrix(rnorm(40, sd = 0.3), nrow = 10), U0)
#' dist_grassmann(X, U0, n = 4, k = 2)
#'
#' @export
dist_grassmann <- function(X, p, n, k) {

    sqrt(rowSums(log_map_grassmann(X, p, n, k)^2))
}

#' Check whether a base point is valid for a Grassmann data set.
#'
#' Checks that \code{p} represents a point of \eqn{Gr(n, k)}, and warns if
#' any row of \code{X} lies on the cut locus of \code{p}, i.e. has a
#' principal angle of \eqn{\pi/2} with it (the smallest singular value of
#' \eqn{U_0^T U} is zero), where the logarithm is not unique.
#'
#' @inheritParams log_map_grassmann
#'
#' @return Logical: \code{TRUE} if no row of \code{X} falls on the cut locus
#'   of \code{p}, \code{FALSE} otherwise. Stops with an error if \code{p}
#'   has the wrong length or non-orthonormal columns.
#'
#' @examples
#' U0 <- diag(4)[, 1:2]
#' X <- exp_map_grassmann(matrix(rnorm(40, sd = 0.3), nrow = 10), U0)
#' is_valid_grassmann(X, U0, n = 4, k = 2)
#'
#' @export
is_valid_grassmann <- function(X, p, n, k) {

    check_grassmann_dims(n, k)
    U0 <- as_stiefel(p, n, k)
    X <- as.matrix(X)
    if (ncol(X) != n * k) {
        stop(sprintf("X must have n * k = %d columns.", n * k))
    }

    cut_locus_idx <- apply(X, 1, function(row) {
        min(svd(crossprod(U0, matrix(row, n, k)), nu = 0, nv = 0)$d) < 1e-10
    })
    n_cut <- sum(cut_locus_idx)
    if (n_cut > 0) {
        warning(sprintf("%d samples detected on the cut locus. Consider modifying the base point.", n_cut))
    }

    invisible(n_cut == 0)
}

#' Find a Grassmann base point via the projection mean.
#'
#' Closed-form, data-centered base point for \eqn{Gr(n, k)}: the subspace
#' spanned by the \code{k} leading eigenvectors of the mean orthogonal
#' projector \eqn{\bar P = N^{-1} \sum_i U_i U_i^T}. This is the extrinsic
#' (chordal) mean of the sample, which is close to the Fréchet mean when
#' the data are concentrated. It is not guaranteed to avoid the cut locus
#' of every sample point, which \code{\link{choose_base_point}} checks
#' afterward.
#'
#' @param X A matrix or data frame of Grassmann data; rows are samples,
#'   with \code{n * k} columns (see \code{\link{manifold_grassmann}}).
#' @inheritParams manifold_grassmann
#'
#' @return Numeric vector of length \code{n * k}: the chosen base point, as
#'   \code{as.vector(U0)}.
#'
#' @examples
#' X <- exp_map_grassmann(matrix(rnorm(40, sd = 0.3), nrow = 10), diag(4)[, 1:2])
#' auto_base_point_grassmann(X, n = 4, k = 2)
#'
#' @export
auto_base_point_grassmann <- function(X, n, k) {

    check_grassmann_dims(n, k)
    X <- as.matrix(X)
    P <- matrix(0, n, n)
    for (i in seq_len(nrow(X))) P <- P + tcrossprod(matrix(X[i, ], n, k))
    as.vector(eigen(P / nrow(X), symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE])
}

#' Find a Grassmann base point via multi-start optimization.
#'
#' Grassmann case of \code{\link{optimize_base_point}}: numerically
#' minimizes the total Grassmann distance (\code{\link{dist_grassmann}})
#' from a candidate point to every row of \code{X}. The search is
#' unconstrained, over the tangent space at the projection mean
#' (\code{\link{auto_base_point_grassmann}}), with restarts drawn from a
#' Gaussian scaled to the spread of the data. See
#' \code{\link{choose_base_point}}.
#'
#' @inheritParams optimize_base_point_torus
#' @param X A matrix or data frame of Grassmann data; rows are samples,
#'   with \code{n * k} columns (see \code{\link{manifold_grassmann}}).
#' @inheritParams manifold_grassmann
#'
#' @return A list with three elements: \code{p} (the best point found
#'   across all restarts, as a vector of length \code{n * k}), \code{value}
#'   (its total Grassmann distance to the rows of \code{X}), and
#'   \code{convergence} (the \code{\link[stats]{optim}} convergence code of
#'   the best restart; \code{0} indicates successful convergence).
#'
#' @examples
#' X <- exp_map_grassmann(matrix(rnorm(40, sd = 0.3), nrow = 10), diag(4)[, 1:2])
#' optimize_base_point_grassmann(X, n = 4, k = 2, n_start = 3)
#'
#' @export
optimize_base_point_grassmann <- function(X, n, k, n_start = 10, maxit = 200, seed = 42) {

    check_grassmann_dims(n, k)
    X <- as.matrix(X)
    d <- (n - k) * k
    center <- matrix(auto_base_point_grassmann(X, n, k), n, k)
    to_point <- function(v) as.vector(exp_map_grassmann(matrix(v, nrow = 1), center))
    spread <- sqrt(mean(dist_grassmann(X, center, n, k)^2) / d)

    # Optimize over the tangent coordinates at the projection mean, which
    # reach every subspace (the exponential map is surjective)
    optimize_base_point(X = X,
                        dist = function(X, v) dist_grassmann(X, to_point(v), n, k),
                        init = function(m) matrix(stats::rnorm(m * d, sd = spread), ncol = d),
                        project = to_point,
                        n_start = n_start,
                        maxit = maxit,
                        seed = seed)
}
