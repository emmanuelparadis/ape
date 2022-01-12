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
  trees[] <- lapply(twoTrees, function (x) x)
  expect_equal(trees, twoTrees)
  
  nexFile <- tempfile()
  write.nexus(twoTrees, file = nexFile)
  # trees read from a nexus file have a $TipLabel attribute
  fromNexus <- read.nexus(nexFile)
  unlink(nexFile) # delete temporary file
  fromNexus[] <- twoTrees
  expect_equal(fromNexus, twoTrees)
})
