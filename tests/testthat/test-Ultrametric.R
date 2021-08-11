test_that("ultrametric trees", {
  set.seed(0)
  N = 100
  n = c(5, 10, 20, 50, 100)
  for (k in n) {
    expect_equal(0, sum(!replicate(N, is.ultrametric(rcoal(k)))))
    expect_equal(0, sum(replicate(N, is.ultrametric(rtree(k)))))
  }
})
