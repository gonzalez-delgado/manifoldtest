test_that("is_in_torus accepts [0,1) and rejects anything else", {
  expect_true(is_in_torus(matrix(c(0, 0.5, 0.999), ncol = 1)))
  expect_true(is_in_torus(data.frame(a = c(0.1, 0.2), b = c(0.3, 0.4))))
  expect_false(is_in_torus(matrix(c(0, 1), ncol = 1)))
  expect_false(is_in_torus(matrix(c(0.5, -0.1), ncol = 1)))
})

test_that("log_map_torus keeps the shape of X and lands in [-0.5, 0.5)", {
  set.seed(1)
  X <- matrix(runif(60), nrow = 20, ncol = 3)
  Z <- log_map_torus(X, p = c(0.1, 0.5, 0.9))

  expect_true(is.matrix(Z))
  expect_equal(dim(Z), dim(X))
  expect_true(all(Z >= -0.5 & Z < 0.5))
})

test_that("log_map_torus sends the cut locus to -0.5", {
  # Exact binary fractions, so the offsets are exactly 0.5
  X <- matrix(c(0.5, 0.75), nrow = 1)
  expect_equal(log_map_torus(X, p = c(0, 0.25)), matrix(-0.5, 1, 2))
})

test_that("log_map_torus gives the shortest signed offset, wrapping around", {
  X <- matrix(c(0.95, 0.05, 0.30), ncol = 1)
  expect_equal(log_map_torus(X, p = 0.05)[, 1], c(-0.1, 0, 0.25))
  expect_equal(log_map_torus(X, p = 0.95)[, 1], c(0, 0.1, 0.35))
})

test_that("log_map_torus at the base point is the origin", {
  p <- c(0.2, 0.7)
  expect_equal(log_map_torus(matrix(p, nrow = 1), p), matrix(0, 1, 2))
})

test_that("log_map_torus accepts data frames", {
  X <- data.frame(a = c(0.1, 0.9), b = c(0.4, 0.6))
  expect_equal(log_map_torus(X, c(0, 0)), log_map_torus(as.matrix(X), c(0, 0)))
})

test_that("dist_torus uses the wraparound distance per coordinate", {
  expect_equal(dist_torus(matrix(0.95, 1, 1), 0.05), 0.1)
  X <- rbind(c(0.9, 0.1), c(0.1, 0.9), c(0.5, 0.5))
  expect_equal(dist_torus(X, c(0, 0)), c(sqrt(0.02), sqrt(0.02), sqrt(0.5)))
})

test_that("dist_torus matches the norm of the log map", {
  set.seed(2)
  X <- matrix(runif(40), nrow = 20, ncol = 2)
  p <- c(0.3, 0.8)
  expect_equal(dist_torus(X, p), sqrt(rowSums(log_map_torus(X, p)^2)))
})

test_that("is_valid_torus detects points on the cut locus", {
  X <- rbind(c(0.1, 0.2), c(0.6, 0.3))
  expect_true(is_valid_torus(X, c(0.3, 0.9)))

  # (0.6 - 0.1) mod 1 = 0.5: first coordinate of row 2 is on the cut locus
  expect_warning(res <- is_valid_torus(X, c(0.1, 0.9)), "1 samples")
  expect_false(res)
})

test_that("is_valid_torus errors on a base point of the wrong length", {
  X <- matrix(runif(6), ncol = 2)
  expect_error(is_valid_torus(X, c(0.1, 0.2, 0.3)), "mismatches")
})

test_that("largest_gap_point_torus returns the midpoint of the largest circular gap", {
  # Gaps: 0.2, 0.05 and 0.75 (wrapping from 0.35 to 0.1)
  expect_equal(largest_gap_point_torus(c(0.1, 0.3, 0.35)), 0.725)
  expect_equal(largest_gap_point_torus(c(0.4, 0.1)), 0.75)
  expect_equal(largest_gap_point_torus(0.2), 0.7)
  expect_equal(largest_gap_point_torus(numeric(0)), 0)
})

test_that("auto_base_point_torus avoids the cut locus of the data", {
  set.seed(3)
  X <- matrix(runif(150), nrow = 50, ncol = 3)
  p <- auto_base_point_torus(X)

  expect_length(p, 3)
  expect_true(all(p >= 0 & p < 1))
  expect_silent(valid <- is_valid_torus(X, p))
  expect_true(valid)
})

test_that("auto_base_point_torus avoids the cut locus of regularly spaced data", {
  X <- matrix(seq(0, 0.9, by = 0.1), ncol = 1)
  expect_silent(valid <- is_valid_torus(X, auto_base_point_torus(X)))
  expect_true(valid)
})

test_that("optimize_base_point_torus finds the center of a cluster across 0/1", {
  X <- matrix(c(0.97, 0.99, 0.01, 0.03), ncol = 1)
  res <- optimize_base_point_torus(X, n_start = 10)

  expect_named(res, c("p", "value", "convergence"))
  expect_true(all(res$p >= 0 & res$p <= 1))
  expect_lt(min(res$p, 1 - res$p), 0.03)
  expect_equal(res$value, sum(dist_torus(X, res$p)))
})

test_that("optimize_base_point_torus is reproducible for a fixed seed", {
  set.seed(4)
  X <- matrix(runif(40), nrow = 20, ncol = 2)
  expect_equal(optimize_base_point_torus(X, n_start = 3, seed = 7),
               optimize_base_point_torus(X, n_start = 3, seed = 7))
})
