lorentz <- function(X, p) -X[, 1] * p[1] + as.vector(X[, -1, drop = FALSE] %*% p[-1])

test_that("exp_map_hyperbolic lands on the hyperboloid", {
  set.seed(20)
  X <- exp_map_hyperbolic(matrix(rnorm(60, sd = 2), nrow = 20, ncol = 3))
  expect_equal(dim(X), c(20, 4))
  expect_true(is_in_hyperbolic(X))
  expect_equal(lorentz(X, X[1, ])[1], -1)
})

test_that("is_in_hyperbolic rejects points off the upper sheet", {
  expect_true(is_in_hyperbolic(matrix(c(1, 0, 0), nrow = 1)))
  expect_false(is_in_hyperbolic(matrix(c(-1, 0, 0), nrow = 1)))
  expect_false(is_in_hyperbolic(matrix(c(2, 0, 0), nrow = 1)))
  expect_false(is_in_hyperbolic(matrix(runif(30), nrow = 10)))
  expect_false(is_in_hyperbolic(matrix(1, nrow = 3, ncol = 1)))
})

test_that("log_map_hyperbolic at the origin matches the closed form of the paper", {
  set.seed(21)
  X <- exp_map_hyperbolic(matrix(rnorm(60), nrow = 20, ncol = 3))
  xr <- X[, -1]
  expected <- acosh(X[, 1]) * xr / sqrt(rowSums(xr^2))
  expect_equal(log_map_hyperbolic(X, c(1, 0, 0, 0)), expected)
})

test_that("log_map_hyperbolic sends the base point to the origin", {
  p <- as.vector(exp_map_hyperbolic(matrix(c(0.4, -1.2), nrow = 1)))
  expect_equal(log_map_hyperbolic(matrix(p, nrow = 1), p), matrix(0, 1, 2))
})

test_that("exp_map_hyperbolic and log_map_hyperbolic are inverse at any base point", {
  set.seed(22)
  p <- as.vector(exp_map_hyperbolic(matrix(c(1.5, -0.7, 0.3), nrow = 1)))
  V <- matrix(rnorm(60, sd = 2), nrow = 20, ncol = 3)
  expect_equal(log_map_hyperbolic(exp_map_hyperbolic(V, p), p), V)
})

test_that("dist_hyperbolic equals arccosh of minus the Lorentzian product", {
  set.seed(23)
  X <- exp_map_hyperbolic(matrix(rnorm(40), nrow = 20, ncol = 2))
  p <- as.vector(exp_map_hyperbolic(matrix(c(1, 1), nrow = 1)))
  expect_equal(dist_hyperbolic(X, p), acosh(-lorentz(X, p)))
})

test_that("dist_hyperbolic is accurate for nearby points", {
  p <- c(1, 0, 0)
  x <- exp_map_hyperbolic(matrix(c(1e-9, 0), nrow = 1))
  expect_equal(dist_hyperbolic(x, p), 1e-9, tolerance = 1e-6)
})

test_that("the log map preserves distances to the base point", {
  set.seed(24)
  X <- exp_map_hyperbolic(matrix(rnorm(20, sd = 2), nrow = 10, ncol = 2))
  p <- c(1, 0, 0)
  expect_equal(sqrt(rowSums(log_map_hyperbolic(X, p)^2)), dist_hyperbolic(X, p))
})

test_that("is_valid_hyperbolic accepts any point on the hyperboloid", {
  X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
  expect_silent(valid <- is_valid_hyperbolic(X, c(1, 0, 0)))
  expect_true(valid)
  expect_error(is_valid_hyperbolic(X, c(1, 0)), "mismatches")
  expect_error(is_valid_hyperbolic(X, c(2, 0, 0)), "hyperboloid")
})

test_that("auto_base_point_hyperbolic is on the hyperboloid and near a concentrated sample", {
  set.seed(25)
  center <- as.vector(exp_map_hyperbolic(matrix(c(2, -1), nrow = 1)))
  X <- exp_map_hyperbolic(matrix(rnorm(100, sd = 0.05), nrow = 50, ncol = 2), center)
  p <- auto_base_point_hyperbolic(X)

  expect_true(is_in_hyperbolic(matrix(p, nrow = 1)))
  expect_lt(dist_hyperbolic(matrix(p, nrow = 1), center), 0.05)
})

test_that("optimize_base_point_hyperbolic improves on the centroid", {
  set.seed(26)
  X <- exp_map_hyperbolic(matrix(rexp(60), nrow = 30, ncol = 2))
  res <- optimize_base_point_hyperbolic(X, n_start = 5)
  auto <- auto_base_point_hyperbolic(X)

  expect_named(res, c("p", "value", "convergence"))
  expect_true(is_in_hyperbolic(matrix(res$p, nrow = 1)))
  expect_equal(res$value, sum(dist_hyperbolic(X, res$p)))
  expect_lte(res$value, sum(dist_hyperbolic(X, auto)) + 1e-8)
})

test_that("choose_base_point works on hyperbolic space by name", {
  X <- exp_map_hyperbolic(matrix(rnorm(20), nrow = 10, ncol = 2))
  res <- choose_base_point(X, manifold = "hyperbolic")
  expect_equal(res$p, auto_base_point_hyperbolic(X))
})

test_that("manifold.dcov.test runs on hyperbolic space", {
  set.seed(27)
  U <- matrix(rnorm(150), nrow = 50, ncol = 3)
  X <- exp_map_hyperbolic(U)
  Y <- exp_map_hyperbolic(U + matrix(rnorm(150, sd = 0.3), nrow = 50))
  res <- manifold.dcov.test(X, Y, manifold = "hyperbolic", R = 99, verbose = FALSE)

  expect_lt(res$dcov_test$p.value, 0.05)
  expect_equal(dim(res$X_mapped), c(50, 3))
  expect_error(manifold.dcov.test(matrix(runif(30), 10), X[1:10, ], manifold = "hyperbolic",
                                  verbose = FALSE),
               "X should lie on the hyperbolic")
})
