U0 <- diag(4)[, 1:2]

# Random tangent coordinates at U0 strictly inside the cut domain
rtangent <- function(N, sd = 0.3) {
  B <- matrix(rnorm(4 * N, sd = sd), nrow = N, ncol = 4)
  s1 <- apply(B, 1, function(r) svd(matrix(r, 2, 2))$d[1])
  B[s1 < 1.5, , drop = FALSE]
}

principal_angles <- function(row, U) acos(pmin(svd(crossprod(U, matrix(row, 4, 2)))$d, 1))

test_that("manifold_grassmann validates its dimensions", {
  expect_s3_class(manifold_grassmann(4, 2), "manifold")
  expect_error(manifold_grassmann(3, 3), "1 <= k <= n - 1")
  expect_error(manifold_grassmann(4, 0), "1 <= k <= n - 1")
  expect_error(manifold_grassmann(4.5, 2), "integers")
})

test_that("the Grassmann manifold cannot be retrieved by name", {
  expect_error(get_manifold("grassmann"), "manifold_grassmann\\(n, k\\)")
})

test_that("exp_map_grassmann lands on Gr(n, k)", {
  set.seed(30)
  X <- exp_map_grassmann(rtangent(20), U0)
  expect_equal(ncol(X), 8)
  expect_true(is_in_grassmann(X, 4, 2))
})

test_that("is_in_grassmann rejects non-orthonormal rows and wrong shapes", {
  expect_true(is_in_grassmann(matrix(as.vector(U0), nrow = 1), 4, 2))
  expect_false(is_in_grassmann(matrix(2 * as.vector(U0), nrow = 1), 4, 2))
  expect_false(is_in_grassmann(matrix(as.vector(U0), nrow = 1), 3, 2))
  expect_false(is_in_grassmann(matrix(runif(80), nrow = 10), 4, 2))
})

test_that("exp_map_grassmann and log_map_grassmann are inverse on the cut domain", {
  set.seed(31)
  B <- rtangent(30, sd = 0.5)
  expect_equal(log_map_grassmann(exp_map_grassmann(B, U0), U0, 4, 2), B)
})

test_that("log_map_grassmann works at a generic base point given as a vector", {
  set.seed(32)
  p <- qr.Q(qr(matrix(rnorm(8), 4, 2)))
  B <- rtangent(10)
  X <- exp_map_grassmann(B, as.vector(p), n = 4, k = 2)
  expect_equal(log_map_grassmann(X, as.vector(p), 4, 2), B)
  expect_equal(log_map_grassmann(matrix(as.vector(p), 1), p, 4, 2), matrix(0, 1, 4))
})

test_that("dist_grassmann is the norm of the principal angles", {
  set.seed(33)
  X <- exp_map_grassmann(rtangent(20, sd = 0.6), U0)
  expected <- apply(X, 1, function(r) sqrt(sum(principal_angles(r, U0)^2)))
  expect_equal(dist_grassmann(X, U0, 4, 2), expected)
})

test_that("dist_grassmann does not depend on the Stiefel representative", {
  set.seed(34)
  X <- exp_map_grassmann(rtangent(10), U0)
  R <- qr.Q(qr(matrix(rnorm(4), 2, 2)))
  XR <- t(apply(X, 1, function(r) as.vector(matrix(r, 4, 2) %*% R)))
  expect_equal(dist_grassmann(XR, U0, 4, 2), dist_grassmann(X, U0, 4, 2))
  expect_equal(log_map_grassmann(XR, U0, 4, 2), log_map_grassmann(X, U0, 4, 2))
})

test_that("log_map_grassmann matches the horizontal lift of Bendokat et al.", {
  # Subspace at principal angles (0.3, 0.7) from span(e1, e2)
  X <- matrix(c(cos(0.3), 0, sin(0.3), 0, 0, cos(0.7), 0, sin(0.7)), nrow = 1)
  coords <- log_map_grassmann(X, U0, 4, 2)
  expect_equal(sort(svd(matrix(coords, 2, 2))$d), c(0.3, 0.7))
})

test_that("log_map_grassmann handles Gr(n, 1) and Gr(n, n - 1)", {
  set.seed(35)
  e1 <- c(1, 0, 0)
  X1 <- exp_map_grassmann(matrix(rnorm(10, sd = 0.3), 5, 2), matrix(e1, 3, 1))
  expect_equal(dim(log_map_grassmann(X1, e1, 3, 1)), c(5, 2))
  p2 <- diag(3)[, 1:2]
  X2 <- exp_map_grassmann(matrix(rnorm(10, sd = 0.3), 5, 2), p2)
  expect_equal(dim(log_map_grassmann(X2, p2, 3, 2)), c(5, 2))
})

test_that("is_valid_grassmann detects subspaces orthogonal to the base point", {
  X <- rbind(as.vector(U0), c(0, 0, 1, 0, 0, 0, 0, 1))
  expect_warning(res <- is_valid_grassmann(X, U0, 4, 2), "1 samples")
  expect_false(res)
  expect_silent(valid <- is_valid_grassmann(X[1, , drop = FALSE], U0, 4, 2))
  expect_true(valid)
  expect_error(is_valid_grassmann(X, 2 * U0, 4, 2), "orthonormal")
  expect_error(is_valid_grassmann(X, U0[, 1], 4, 2), "mismatches")
})

test_that("auto_base_point_grassmann is near a concentrated sample", {
  set.seed(36)
  center <- qr.Q(qr(matrix(rnorm(8), 4, 2)))
  X <- exp_map_grassmann(matrix(rnorm(200, sd = 0.03), 50, 4), center)
  p <- auto_base_point_grassmann(X, 4, 2)

  expect_true(is_in_grassmann(matrix(p, nrow = 1), 4, 2))
  expect_lt(dist_grassmann(matrix(as.vector(center), 1), p, 4, 2), 0.05)
})

test_that("optimize_base_point_grassmann improves on the projection mean", {
  set.seed(37)
  X <- exp_map_grassmann(rtangent(30, sd = 0.5), U0)
  res <- optimize_base_point_grassmann(X, 4, 2, n_start = 3)
  auto <- auto_base_point_grassmann(X, 4, 2)

  expect_true(is_in_grassmann(matrix(res$p, nrow = 1), 4, 2))
  expect_equal(res$value, sum(dist_grassmann(X, res$p, 4, 2)))
  expect_lte(res$value, sum(dist_grassmann(X, auto, 4, 2)) + 1e-8)
})

test_that("manifold.dcov.test runs on the Grassmann manifold", {
  set.seed(38)
  B <- rtangent(60, sd = 0.4)
  X <- exp_map_grassmann(B, U0)
  Y <- exp_map_grassmann(B + matrix(rnorm(length(B), sd = 0.1), nrow(B)), U0)
  m <- manifold_grassmann(4, 2)
  res <- manifold.dcov.test(X, Y, manifold = m, center_outward = TRUE, R = 99, verbose = FALSE)

  expect_lt(res$dcov_test$p.value, 0.05)
  expect_equal(ncol(res$X_mapped), 4)
  expect_true(all(sqrt(rowSums(res$X_mapped^2)) < 1))
  expect_error(manifold.dcov.test(X + 1, Y, manifold = m, verbose = FALSE), "X should lie on")
})
