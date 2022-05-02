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
  fromNexus <- read.nexus(nexFile)
  unlink(nexFile) # delete temporary file
  
  # trees read from a nexus file have a $TipLabel attribute
  tipLabel <- attr(fromNexus, "TipLabel")
  expect_equal(length(tipLabel), 6)
  attr(fromNexus, "names") <- NULL
  withoutTL <- c(fromNexus[[1]], fromNexus[[2]])
  bothLabels <- structure(withoutTL, TipLabel = tipLabel)
  
  trees <- fromNexus
  trees[] <- fromNexus
  expect_equal(trees, fromNexus)
  
  trees[] <- fromNexus[2:1]
  expect_equal(trees, rev(fromNexus))
  
  trees[1] <- trees[2]
  expect_equal(trees, fromNexus[c(1, 1)])
  Capitalize <- function (tree) {
    tree$tip.label <- toupper(tree$tip.label)
    # Return:
    tree
  }
  trees <- fromNexus
  trees[] <- lapply(trees, Capitalize)
  expect_equal(trees,
               structure(lapply(fromNexus, Capitalize), class = 'multiPhylo'))
})
