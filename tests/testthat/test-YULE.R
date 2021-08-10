test_that("Random Yule trees", {
  N = 1000
  lambda = 0.05
  Tmax = 50
  threshold = c(0.8, 1.2)
  x <- as.integer(replicate(floor(N/2), balance(rlineage(lambda, 0, Tmax))[1, ]))
  mx <- max(x)
  O <- tabulate(x, mx)
  P <- length(x) * dyule(1:mx, lambda, Tmax)
  r <- cor(P, O)
  expect_gt(r, threshold[1])
  expect_lt(r, threshold[2])
})
