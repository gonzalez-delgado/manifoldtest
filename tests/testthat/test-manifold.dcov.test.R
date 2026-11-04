test_that("manifold.dcov.test returns a p-value in [0, 1] and mapped data", {
  set.seed(12)
  X <- matrix(runif(90), nrow = 30, ncol = 3)
  Y <- matrix(runif(60), nrow = 30, ncol = 2)
  res <- manifold.dcov.test(X, Y, R = 99, center_outward = TRUE, verbose = FALSE)

  expect_named(res, c("dcov_test", "px_base", "py_base", "X_mapped", "Y_mapped"))
  expect_s3_class(res$dcov_test, "htest")
  expect_gte(res$dcov_test$p.value, 0)
  expect_lte(res$dcov_test$p.value, 1)
  expect_equal(dim(res$X_mapped), c(30, 3))
  expect_equal(dim(res$Y_mapped), c(30, 2))
  expect_true(all(sqrt(rowSums(res$X_mapped^2)) < 1))
})

test_that("manifold.dcov.test runs DistCov by default", {
  set.seed(19)
  X <- matrix(runif(40), nrow = 20, ncol = 2)
  Y <- matrix(runif(40), nrow = 20, ncol = 2)
  res <- manifold.dcov.test(X, Y, R = 19, verbose = FALSE)
  expect_equal(res$X_mapped, log_map_torus(X, res$px_base))
})

test_that("manifold.dcov.test without center-outward uses the raw log map", {
  set.seed(13)
  X <- matrix(runif(40), nrow = 20, ncol = 2)
  Y <- matrix(runif(40), nrow = 20, ncol = 2)
  res <- manifold.dcov.test(X, Y, R = 99, center_outward = FALSE, verbose = FALSE)

  expect_equal(res$X_mapped, log_map_torus(X, res$px_base))
  expect_equal(res$Y_mapped, log_map_torus(Y, res$py_base))
})

test_that("manifold.dcov.test keeps a valid user-supplied base point", {
  X <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.15, 0.25, 0.35, 0.45), ncol = 2)
  Y <- matrix(c(0.6, 0.7, 0.8, 0.9), ncol = 1)
  res <- manifold.dcov.test(X, Y, px_base = c(0.2, 0.3), py_base = 0.75,
                            R = 19, verbose = FALSE)
  expect_equal(res$px_base, c(0.2, 0.3))
  expect_equal(res$py_base, 0.75)
})

test_that("manifold.dcov.test rejects data outside [0,1)", {
  X <- matrix(runif(20), ncol = 2)
  expect_error(manifold.dcov.test(X + 1, X, verbose = FALSE), "X should lie on the torus")
  expect_error(manifold.dcov.test(X, X - 1, verbose = FALSE), "Y should lie on the torus")
})

test_that("manifold.dcov.test aligns unequal sample sizes", {
  set.seed(14)
  X <- matrix(runif(60), nrow = 30, ncol = 2)
  Y <- matrix(runif(40), nrow = 20, ncol = 2)
  expect_message(res <- manifold.dcov.test(X, Y, R = 19, verbose = FALSE),
                 "aligned")
  expect_equal(nrow(res$X_mapped), 20)
  expect_equal(nrow(res$Y_mapped), 20)
})

test_that("manifold.dcov.test detects strong dependence across the 0/1 seam", {
  set.seed(15)
  X <- matrix(runif(100), nrow = 50, ncol = 2)
  # Y is a noisy copy of X, wrapped back onto the torus
  Y <- (X + matrix(rnorm(100, sd = 0.02), nrow = 50)) %% 1
  res <- manifold.dcov.test(X, Y, R = 199, verbose = FALSE)
  expect_lt(res$dcov_test$p.value, 0.01)
})

test_that("manifold.dcov.test p-values are roughly uniform under independence", {
  skip_on_cran()

  n_rep <- 100
  set.seed(16)
  data <- replicate(n_rep, list(X = matrix(runif(60), 30, 2),
                                Y = matrix(runif(60), 30, 2)),
                    simplify = FALSE)
  pvals <- vapply(data, function(d)
    manifold.dcov.test(d$X, d$Y, R = 99, verbose = FALSE)$dcov_test$p.value,
    numeric(1))

  # Under H0, rejections at 5% ~ Binomial(100, 0.05): P(> 12) < 0.001
  expect_lte(sum(pvals < 0.05), 12)
  expect_gt(mean(pvals), 0.35)
  expect_lt(mean(pvals), 0.65)
})

test_that("manifold.dcov.test accepts vectors as one-dimensional samples", {
  set.seed(18)
  x <- runif(30)
  y <- (x + rnorm(30, sd = 0.05)) %% 1
  res <- manifold.dcov.test(x, y, R = 99, verbose = FALSE)
  expect_equal(dim(res$X_mapped), c(30, 1))
  expect_lt(res$dcov_test$p.value, 0.05)
})
