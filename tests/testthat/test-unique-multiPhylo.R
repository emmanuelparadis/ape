test_that("unique.multiPhylo()", {
  tree <- read.tree(text = '(a, b);')
  expect_equal(unique(c(tree, tree)), structure(c(tree), old.index = rep(1L, 2)))
  expect_equal(unique(c(tree)), structure(c(tree), old.index = rep(1L, 1)))
  expect_equal(unique(c(tree)[-1]), structure(c(tree)[-1]))
})
