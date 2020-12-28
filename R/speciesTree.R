## speciesTree.R (2013-08-12)

##   Species Trees

## Copyright 2010-2013 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

speciesTree <- function(x, FUN = min)
### FUN = min => MAXTREE (Liu et al. 2010)
### FUN = sum => shallowest divergence (Maddison & Knowles 2006)
{
    test.ultra <- which(!unlist(lapply(x, is.ultrametric)))
    if (length(test.ultra))
        stop(paste("the following trees were not ultrametric:\n",
                   paste(test.ultra, collapse = " ")))

    Ntree <- length(x)
    D <- lapply(x, cophenetic.phylo)
    nms <- rownames(D[[1]])
    n <- length(nms)
    M <- matrix(0, n*(n - 1)/2, Ntree)
    for (i in 1:Ntree) M[, i] <- as.dist(D[[i]][nms, nms])
    Y <- apply(M, 1, FUN)
    attributes(Y) <- list(Size = n, Labels = nms, Diag = FALSE,
                          Upper = FALSE, class = "dist")
    as.phylo(hclust(Y, "single"))
}
