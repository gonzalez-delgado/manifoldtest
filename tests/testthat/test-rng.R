# The package uses fixed seeds internally for reproducibility; these must not
# reset the caller's random number generator.

expect_rng_untouched <- function(expr) {
  set.seed(123)
  before <- .Random.seed
  force(expr)
  expect_identical(.Random.seed, before)
}

test_that("build_ball_grid does not alter the global RNG", {
  expect_rng_untouched(build_ball_grid(n = 10, d = 3, nR = 3, nS = 3, n0 = 1))
})

test_that("align_sample_sizes does not alter the global RNG", {
  X <- matrix(1:20, nrow = 10)
  Y <- matrix(1:12, nrow = 6)
  expect_rng_untouched(suppressMessages(align_sample_sizes(X, Y)))
})

test_that("optimize_base_point_torus does not alter the global RNG", {
  X <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 2)
  expect_rng_untouched(optimize_base_point_torus(X, n_start = 2))
})

test_that("repeated tests on fresh data give different p-values", {
  set.seed(17)
  p <- replicate(5, {
    X <- matrix(runif(60), 20, 3)
    Y <- matrix(runif(60), 20, 3)
    manifold.dcov.test(X, Y, R = 99, verbose = FALSE)$dcov_test$p.value
  })
  expect_gt(length(unique(p)), 1)
})
