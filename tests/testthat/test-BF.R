test_that("Base frequencies calculated", {
  fas1 <- system.file("extdata/input/FASTA/seq1_DNA.fas", package = "ape")
  dna1 <- read.dna(fas1, format = "f")
  BF1 <- base.freq(dna1, TRUE, TRUE)
  out1 <- system.file("extdata/output/BF1.txt", package = "ape")
  BF1.0 <- as.double(read.table(out1, header = TRUE))
  expect_equal(unname(BF1[c("a", "c", "g", "t", "n")]), BF1.0)
})
