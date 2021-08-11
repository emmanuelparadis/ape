test_that("Phylogenetically independent contrasts values agree", {
  treefile <- system.file("extdata/input/Newick/tree_primates.tre", package = "ape")
  datfile <- system.file("extdata/input/Table/data_primates.txt", package = "ape")
  tree.primates <- read.tree(treefile)
  DATA <- read.table(datfile, header = TRUE)
  pic.body <- pic(DATA$body, tree.primates)
  pic.brain <- pic(DATA$brain, tree.primates)
  outfile <- system.file("extdata/output/PIC_primates.txt", package = "ape")
  PIC.0 <- read.table(outfile, header = TRUE)
  ## only 6 digits in PHYLIP's output
  expect_equal(unname(sort(pic.body)), sort(PIC.0$body), tolerance = 1e-5)
  expect_equal(unname(sort(pic.brain)), sort(PIC.0$brain), tolerance = 1e-5)
})
