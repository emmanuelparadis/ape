test_that("write.tree() writes lists", {
  trees <- c('((a,b),(c,d));', '(a,(b,d,c));')
  phy <- read.tree(text = trees)

  expect_equal(write.tree(phy), trees)
  expect_equal(write.tree(unclass(phy)), trees)
  expect_error(write.tree(c('not-a-tree', unclass(phy))))
})
