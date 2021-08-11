test_that("Split frequencies correct", {
  fl <- system.file("extdata/input/Newick/three_unrooted_trees_4tips.tre",
                    package = "ape")
  TR <- read.tree(fl)
  a <- summary(prop.part(TR))[-1]
  b <- bitsplits(TR)$freq
  expect_equal(a, rep(1, 3))
  expect_equal(b, rep(1, 3))
})
