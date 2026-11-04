test_that("get_manifold resolves names and passes manifold objects through", {
  m <- get_manifold("torus")
  expect_s3_class(m, "manifold")
  expect_equal(m$name, "torus")
  expect_identical(get_manifold(m), m)
})

test_that("get_manifold retrieves hyperbolic space by name", {
  expect_equal(get_manifold("hyperbolic")$name, "hyperbolic")
})

test_that("get_manifold rejects unknown or malformed manifolds", {
  expect_error(get_manifold("sphere"), "Unknown manifold 'sphere'")
  expect_error(get_manifold(c("torus", "torus")), "single character string")
  expect_error(get_manifold(1), "single character string")
})

test_that("manifold_torus exposes the torus functions and base point methods", {
  m <- manifold_torus()
  expect_identical(m$log_map, log_map_torus)
  expect_identical(m$dist, dist_torus)
  expect_identical(m$is_valid, is_valid_torus)
  expect_identical(m$contains, is_in_torus)
  expect_equal(names(m$base_point), c("auto", "optimize"))
})

test_that("new_manifold validates its arguments", {
  f <- function(...) NULL
  expect_s3_class(
    new_manifold("m", "R^d", f, f, f, f, list(default = f)),
    "manifold"
  )
  expect_error(new_manifold("m", "R^d", f, f, f, f, list()))
  expect_error(new_manifold("m", "R^d", f, f, f, f, list(f)))
  expect_error(new_manifold("m", "R^d", f, "not a function", f, f, list(a = f)))
})

test_that("print.manifold shows the name and base point methods", {
  expect_output(print(manifold_torus()), "<manifold: torus>")
  expect_output(print(manifold_torus()), "auto, optimize")
})
