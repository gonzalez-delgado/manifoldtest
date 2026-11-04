test_that("choose_base_point uses the 'auto' method by default", {
  set.seed(10)
  X <- matrix(runif(40), nrow = 20, ncol = 2)
  res <- choose_base_point(X)

  expect_named(res, c("p", "value", "convergence"))
  expect_equal(res$p, auto_base_point_torus(X))
  expect_equal(res$value, sum(dist_torus(X, res$p)))
  expect_true(is.na(res$convergence))
})

test_that("choose_base_point 'optimize' is at least as central as 'auto'", {
  set.seed(11)
  X <- matrix(runif(40), nrow = 20, ncol = 2)
  auto <- choose_base_point(X, method = "auto")
  opt  <- choose_base_point(X, method = "optimize", n_start = 5)

  expect_lte(opt$value, auto$value + 1e-8)
})

test_that("choose_base_point rejects unknown methods and empty data", {
  X <- matrix(runif(10), ncol = 2)
  expect_error(choose_base_point(X, method = "foo"))
  expect_error(choose_base_point(matrix(numeric(0), nrow = 5, ncol = 0)),
               "at least 1 column")
})

test_that("choose_base_point works with a custom manifold", {
  m <- new_manifold(
    name = "toy", representation = "R^d",
    contains = function(X) TRUE,
    log_map = function(X, p) sweep(as.matrix(X), 2, p),
    dist = function(X, p) sqrt(rowSums(sweep(as.matrix(X), 2, p)^2)),
    is_valid = function(X, p) TRUE,
    base_point = list(mean = function(X, ...) list(p = colMeans(X), convergence = NA))
  )
  X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  expect_equal(choose_base_point(X, manifold = m)$p, c(2, 5))
})
