test_that("Type I error rate of the Mantel test", {
  set.seed(0)
  N = 100 * 200
  n = 10
  falsePositiveRate = 0.01

  rmat <- function(n) {
    x <- runif(n * (n - 1) / 2)
    m <- matrix(0, n, n)
    m[lower.tri(m)] <- x
    m <- t(m)
    m[lower.tri(m)] <- x
    m
  }
  res <- numeric()
  res <- replicate(N, {
    ma <- rmat(n)
    mb <- rmat(n)
    mantel.test(ma, mb)$p
  })
  expect_false(anyNA(res))
  expect_lte(pbinom(sum(res < 0.05), length(res), 0.05), 1 - falsePositiveRate)
})
