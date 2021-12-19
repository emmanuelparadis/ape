test_that("[<-.multiPhylo works", {
  set.seed(0)
  twoTrees <- c(rtree(6), rtree(6))
  trees <- twoTrees
  trees[] <- trees
  expect_equal(trees, twoTrees)
  trees[] <- trees[2:1]
  expect_equal(trees, rev(twoTrees))
  trees[1] <- trees[2]
  expect_equal(trees, c(twoTrees[1], twoTrees[1]))
})
