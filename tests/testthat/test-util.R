test_that("align_sample_sizes leaves equal sizes unchanged", {
  X <- matrix(1:6, nrow = 3)
  Y <- matrix(7:12, nrow = 3)
  expect_silent(res <- align_sample_sizes(X, Y))
  expect_identical(res, list(X = X, Y = Y))
})

test_that("align_sample_sizes subsamples to the smaller size", {
  X <- matrix(1:20, nrow = 10)
  Y <- matrix(1:12, nrow = 6)
  expect_message(res <- align_sample_sizes(X, Y), "inconsistent")

  expect_equal(nrow(res$X), 6)
  expect_equal(nrow(res$Y), 6)
  expect_equal(ncol(res$X), 2)
  # Rows of X are a subset of the original rows, without repetition
  expect_true(all(res$X[, 1] %in% X[, 1]))
  expect_false(anyDuplicated(res$X[, 1]) > 0)
})

test_that("align_sample_sizes is reproducible for a fixed seed", {
  X <- matrix(1:20, nrow = 10)
  Y <- matrix(1:12, nrow = 6)
  expect_equal(suppressMessages(align_sample_sizes(X, Y, seed = 1)),
               suppressMessages(align_sample_sizes(X, Y, seed = 1)))
})

test_that("optimize_base_point checks the number of restart points", {
  X <- matrix(runif(10), ncol = 2)
  expect_error(
    optimize_base_point(X, dist = dist_torus,
                        init = function(n) matrix(0.5, nrow = 1, ncol = 2),
                        n_start = 3),
    "returned 1 restart points instead of 3"
  )
})
