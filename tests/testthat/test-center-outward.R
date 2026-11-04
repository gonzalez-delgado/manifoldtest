test_that("build_ball_grid has the requested size and radial structure", {
  grid <- build_ball_grid(n = 14, d = 2, nR = 3, nS = 4, n0 = 2)
  radii <- sqrt(rowSums(grid^2))

  expect_equal(dim(grid), c(14, 2))
  expect_equal(sum(radii == 0), 2)
  expect_equal(sort(unique(round(radii, 12))), c(0, 0.25, 0.5, 0.75))
  expect_equal(as.vector(table(round(radii[radii > 0], 12))), c(4, 4, 4))
})

test_that("build_ball_grid puts every point of a shell at that shell's radius", {
  for (d in 2:4) {
    grid <- build_ball_grid(n = 16, d = d, nR = 3, nS = 5, n0 = 1)
    radii <- sqrt(rowSums(grid^2))
    expect_equal(radii, c(0, rep(c(0.25, 0.5, 0.75), each = 5)))
  }
})

test_that("build_ball_grid stays strictly inside the unit ball for d = 1, 2, 3", {
  for (d in 1:3) {
    grid <- build_ball_grid(n = 13, d = d, nR = 3, nS = 4, n0 = 1)
    expect_equal(dim(grid), c(13, d))
    expect_true(all(sqrt(rowSums(grid^2)) < 1))
  }
})

test_that("build_ball_grid is reproducible for d > 2", {
  expect_equal(build_ball_grid(10, 3, 3, 3, 1, seed = 1),
               build_ball_grid(10, 3, 3, 3, 1, seed = 1))
})

test_that("optimal_assignment returns a permutation", {
  set.seed(5)
  Z <- matrix(rnorm(40), nrow = 20, ncol = 2)
  grid <- build_ball_grid(n = 20, d = 2, nR = 4, nS = 5, n0 = 0)
  a <- optimal_assignment(Z, grid)

  expect_length(a, 20)
  expect_setequal(a, 1:20)
})

test_that("optimal_assignment minimizes the total squared cost", {
  set.seed(6)
  Z <- matrix(rnorm(10), nrow = 5, ncol = 2)
  grid <- build_ball_grid(n = 5, d = 2, nR = 2, nS = 2, n0 = 1)
  cost <- function(perm) sum(rowSums((Z - grid[perm, , drop = FALSE])^2))

  # Brute force over all 5! permutations
  perms <- function(v) {
    if (length(v) <= 1) return(list(v))
    do.call(c, lapply(seq_along(v), function(i)
      lapply(perms(v[-i]), function(r) c(v[i], r))))
  }
  best <- min(vapply(perms(1:5), cost, numeric(1)))

  expect_equal(cost(optimal_assignment(Z, grid)), best)
})

test_that("center_outward_transform maps into the unit ball", {
  set.seed(7)
  for (d in 1:3) {
    Z <- matrix(rnorm(50 * d, sd = 100), nrow = 50, ncol = d)
    F <- center_outward_transform(Z)
    expect_equal(dim(F), c(50, d))
    expect_true(all(sqrt(rowSums(F^2)) < 1))
  }
})

test_that("center_outward_transform uses every grid point exactly once", {
  set.seed(8)
  Z <- matrix(rnorm(40), nrow = 20, ncol = 2)
  F <- center_outward_transform(Z)
  # n = 20 -> nR = 4, nS = 5, n0 = 0
  grid <- build_ball_grid(n = 20, d = 2, nR = 4, nS = 5, n0 = 0)

  key <- function(M) sort(apply(round(M, 10), 1, paste, collapse = ","))
  expect_equal(key(F), key(grid))
})

test_that("center_outward_transform sends outlying points to the outer shell", {
  set.seed(9)
  Z <- rbind(matrix(rnorm(48, sd = 0.1), ncol = 2), c(50, 50))
  F <- center_outward_transform(Z)
  radii <- sqrt(rowSums(F^2))
  expect_equal(radii[25], max(radii))
})
