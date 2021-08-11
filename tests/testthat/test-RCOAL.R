test_that("Random coalescent trees", {
  set.seed(0)
  bound_p <- 0.995
  BOUND <- qnorm(bound_p)
  N <- 100
  tree.sizes <- c(5, 10, 20, 50, 75, 100)
  failures <- 0L
  reps <- 200L
  for (i in seq_len(reps)) {
    for (n in tree.sizes) {
      k <- 2:n
      expected.mean <- 2 * sum(1/(k * (k - 1)))
      expected.var <- 4 * sum(1/(k * (k - 1))^2)
      x <- replicate(N, branching.times(rcoal(n))[1])
      expect_false(anyNA(x))
      failures <- failures + sum(abs((mean(x) - expected.mean) * sqrt(N/expected.var)) > BOUND)
    }
  }
  n <- reps * length(tree.sizes) * N
  q <- bound_p
  p <- 1 - q
  expected_fails <- n * p
  fails_var <- n * p * q
  fails_sd <- sqrt(fails_var)
  # if working: 99.7% chance of lying within 3sd of mean
  fails_bound <- 3 * fails_sd
  expect_lt(failures, expected_fails + fails_bound)
  #expect_gt(failures, expected_fails - fails_bound)
})
