test_that("reorder.phylo works", {
  set.seed(0)
  Nmin = 3L
  Nmax = 1000L
  ProbRooted = 0.5
  ProbMultichotomy = 0.5
  nrep = 1e4L

  Ntip <- Nnode <- integer(nrep)
  Rooted <- Test1 <- Test2 <- logical(nrep)

  pm <- runif(nrep) > ProbMultichotomy

  for (i in 1:nrep) {
    n <- sample(Nmin:Nmax, 1L)
    rooted <- sample(c(TRUE, FALSE), 1L, prob = c(ProbRooted, 1 - ProbRooted))
    tr <- rtree(n, rooted)
    if (pm[i]) {
      if (n == 3 && !rooted) break
      INTS <- which(tr$edge[, 2L] > n)
      m <- length(INTS)
      if (!m) break
      k <- sample(INTS, ceiling(ProbMultichotomy * m))
      tr$edge.length[k] <- 0
      tr <- di2multi(tr)
    }
    Ntip[i] <- Ntip(tr)
    Nnode[i] <- Nnode(tr)
    Rooted[i] <- rooted
    expect_identical(reorder(reorder(tr, "pr"))$edge, tr$edge)
    expect_identical(reorder(reorder(tr, "po"))$edge, tr$edge)
  }
})
