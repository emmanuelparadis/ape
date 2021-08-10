test_that("Variance-covariance under Brownian motion", {
  tr <- compute.brtime(stree(5, "l"), 4:1)
  vcvape <- vcv(tr)
  expected.vcv <- diag(4, 5, 5)
  expected.vcv[lower.tri(expected.vcv)] <- offdiag <- rep(0:3, 4:1)
  expected.vcv <- t(expected.vcv)
  expected.vcv[lower.tri(expected.vcv)] <- offdiag
  dimnames(expected.vcv) <- list(paste0('t', 1:5), paste0('t', 1:5))
  expect_equal(vcvape, expected.vcv)
})
