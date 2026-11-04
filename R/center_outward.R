#' Build a regular grid on the unit ball
#'
#' Following Hallin et al. (2021) and the OT2.R reference code, builds a
#' grid of \code{n} points on the unit ball \eqn{S^d} used as the target of
#' the center-outward transform (\code{\link{center_outward_transform}}).
#' The grid consists of:
#' \itemize{
#'   \item \code{n0} copies of the origin;
#'   \item \code{nR} radial levels \eqn{r_k = k/(n_R+1)}, \eqn{k = 1,
#'     \ldots, n_R}, each carrying the same \code{nS} angular directions;
#' }
#' so that \code{n0 + nR * nS = n}. The \code{nS} directions are
#' quasi-uniform unit vectors on \eqn{S^{d-1}} (equally spaced angles for
#' \code{d = 2}, normalized Gaussian draws for \code{d > 2}, as in OT2.R).
#'
#' @param n Integer: total number of grid points.
#' @param d Integer: dimension of the ball.
#' @param nR Integer: number of radial levels.
#' @param nS Integer: number of angular directions per radial level.
#' @param n0 Integer: number of copies of the origin, so that
#'   \code{n0 + nR * nS = n}.
#' @param seed Integer random seed for the quasi-uniform directions when
#'   \code{d > 2} (unused for \code{d <= 2}, which are deterministic). The
#'   grid is reproducible for a given \code{(nS, seed)}, since \code{nS} is
#'   folded into the seed actually used. Defaults to \code{42}.
#'
#' @return Numeric matrix (\code{n} x \code{d}): rows are the grid points
#'   in the unit ball \eqn{S^d}.
#'
#' @examples
#' # n = n0 + nR * nS = 2 + 3 * 4
#' grid <- build_ball_grid(n = 14, d = 2, nR = 3, nS = 4, n0 = 2)
#' dim(grid)
#'
#' @export
build_ball_grid <- function(n, d, nR, nS, n0, seed = 42) {
  
  # Directions: nS unit vectors on S^{d-1}
  U <- matrix(NA, d, nS)

  if (d == 1) {
    # S^0 = {-1, 1}: alternate the two signs
    U[1, ] <- ifelse((1:nS) %% 2 == 1, 1, -1)
  } else if (d == 2) {
    for (j in 1:nS) {
      U[1, j] <- cos(2 * pi * (j - 1) / nS)
      U[2, j] <- sin(2 * pi * (j - 1) / nS)
    }
  } else {
    # Quasi-uniform directions via normalised Gaussian draws (OT2.R);
    # fixed seed makes the grid reproducible for a given (nS, d)
    set.seed(seed + nS)
    UU <- matrix(stats::rnorm(nS * d), nS, d)
    UU <- UU / sqrt(rowSums(UU^2))
    U  <- t(UU)
  }

  # Radial lengths: n0 zeros, then radii 1/(nR+1)..nR/(nR+1), each x nS
  length_Grid <- c(rep(0, n0),
                   rep(seq(1 / (nR + 1), nR / (nR + 1), length.out = nR),
                       each = nS))

  # Directions in the same column ordering
  direction_Grid <- cbind(matrix(0, d, n0),
                          matrix(rep(U, nR), nrow = d))

  # Grid points = length * direction  (d x n), transpose to n x d
  return(t(length_Grid * direction_Grid))
}

#' Solve the optimal L2 assignment between data and grid
#'
#' Computes the bijection \eqn{\pi} between the \code{n} sample points and
#' the \code{n} grid points that minimizes
#' \eqn{\sum_i \|Z_i - grid_{\pi(i)}\|^2}. This is solved as an optimal
#' transport problem (network-flow solver from the \code{transport}
#' package), which is substantially faster than the Hungarian algorithm
#' implemented in \code{clue::solve_LSAP}. Used by
#' \code{\link{center_outward_transform}} to assign each observation its
#' rank on the grid built by \code{\link{build_ball_grid}}.
#'
#' @param Z Numeric matrix (\code{n} x \code{d}): Euclidean data.
#' @param grid Numeric matrix (\code{n} x \code{d}): target grid, e.g. from
#'   \code{\link{build_ball_grid}}.
#'
#' @return Integer vector of length \code{n}: the grid row assigned to
#'   each row of \code{Z}.
#'
#' @examples
#' Z <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' grid <- build_ball_grid(n = 10, d = 2, nR = 3, nS = 3, n0 = 1)
#' optimal_assignment(Z, grid)
#'
#' @export
optimal_assignment <- function(Z, grid) {

  n <- nrow(Z)
  # Squared Euclidean cost between all data-grid pairs
  distMat <- as.matrix(stats::dist(rbind(Z, grid)))
  cost    <- distMat[1:n, (n + 1):(2 * n)]^2        
  transport::transport(a=rep(1, n),
                      b=rep(1, n),
                      costm = cost,
                      method = "networkflow")$to
}

#' Center-outward transform
#'
#' Maps \code{d}-dimensional Euclidean data to the open unit ball
#' \eqn{S^d}, following Hallin et al. (2021), by returning the empirical
#' center-outward distribution function
#' \eqn{F^{(n)}: \mathbb{R}^d \to S^d}: an \code{n} x \code{d} matrix whose
#' rows are the grid values (\code{\link{build_ball_grid}}) assigned to
#' each observation by the optimal \eqn{L^2} coupling
#' (\code{\link{optimal_assignment}}). \eqn{F} is bounded in the unit ball,
#' so any test statistic computed on \eqn{F} (e.g. a distance covariance
#' test) requires no moment condition on the original data.
#'
#' @param Z Numeric matrix (\code{n} x \code{d}): Euclidean data, typically
#'   the output of a manifold log map (e.g. \code{\link{log_map_torus}}).
#' @param seed Integer random seed passed to \code{\link{build_ball_grid}}
#'   for the grid's quasi-uniform directions (only used when \code{d > 2}).
#'   Defaults to \code{42}.
#'
#' @return Numeric matrix \code{F} (\code{n} x \code{d}): rows inside the
#'   unit ball \eqn{S^d}, the multidimensional rank of each row of
#'   \code{Z}.
#'
#' @examples
#' Z <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' F <- center_outward_transform(Z)
#' dim(F)
#' max(sqrt(rowSums(F^2))) <= 1
#'
#' @export
center_outward_transform <- function(Z, seed = 42) {

  n <- nrow(Z)
  d <- ncol(Z)

  # Factorisation n = n_R * n_S + n_0 (n_R, n_S both grow with n)
  nR <- max(1, floor(sqrt(n)))
  nS <- floor(n / nR)
  n0 <- n - nR * nS
  if (nS < 1) { nR <- n; nS <- 1; n0 <- 0 }

  # Regular grid on the unit ball
  grid <- build_ball_grid(n=n,
                          d=d,
                          nR=nR,
                          nS=nS,
                          n0=n0,
                          seed = seed)

  # Optimal L2 assignment (fast network-flow solver)
  assignment <- optimal_assignment(Z=Z,
                                  grid=grid)

  # Empirical center-outward distribution function
  return(grid[assignment, , drop = FALSE])                  
}
