test_that("Branching times calculated", {
  nwk1 <- system.file("extdata/input/Newick/tree1_Newick.tre",
                      package = "ape")
  tr1 <- read.tree(nwk1)
  bt1 <- branching.times(tr1)
  bt1.0 <- scan(system.file("extdata/output/bt1.txt", package = "ape"),
                sep = "\n", quiet = TRUE)
  expect_equal(unname(bt1), bt1.0)
})
