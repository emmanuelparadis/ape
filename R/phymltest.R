## phymltest.R (2014-11-07)

##   Fits a Bunch of Models with PhyML

## Copyright 2004-2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.phymltest.model <-
    c("JC69", "JC69+I", "JC69+G", "JC69+I+G",
      "K80", "K80+I", "K80+G", "K80+I+G",
      "F81", "F81+I", "F81+G", "F81+I+G",
      "F84", "F84+I", "F84+G", "F84+I+G",
      "HKY85", "HKY85+I", "HKY85+G", "HKY85+I+G",
      "TN93", "TN93+I", "TN93+G", "TN93+I+G",
      "GTR", "GTR+I", "GTR+G", "GTR+I+G")

.phymltest.nfp <-
    c(1, 2, 2, 3, 2, 3, 3, 4, 4, 5, 5, 6, 5, 6, 6, 7,
      5, 6, 6, 7, 6, 7, 7, 8, 9, 10, 10, 11)

phymltest <- function(seqfile, format = "interleaved", itree = NULL,
                      exclude = NULL, execname = NULL, append = TRUE)
{
    os <- Sys.info()[1]
    ## default names of PhyML:
    if (is.null(execname)) {
        execname <- switch(os, "Linux" = {
            ## PhyML location for Debian and Fedora packages and maybe for other distributions (fix by Dylan A\"issi)
            if (file.exists("/usr/bin/phyml")) "/usr/bin/phyml" else "phyml_3.0.1_linux32"
	}, "Darwin" = "phyml_3.0.1_macintel", "Windows" = "phyml_3.0.1_win32")
    }
    if (is.null(execname))
        stop("you must give an executable file name for PHYML")
    N <- length(.phymltest.model)
    format <- match.arg(format, c("interleaved", "sequential"))
    fmt <- rep("", N)
    if (format != "interleaved") fmt[] <- "-q"
    boot <- rep("-b 0", N) # to avoid any testing
    mdl <- paste("-m", rep(c("JC69", "K80", "F81", "F84", "HKY85", "TN93", "GTR"), each = 4)) # fix by Luiz Max Fagundes de Carvalho
    tstv <- rep("-t e", N) # ignored by PhyML with JC69 or F81
    inv <- rep(c("", "-v e"), length.out = N)
    ## no need to use the -c option of PhyML (4 categories by default if '-a e' is set):
    alpha <- rep(rep(c("-c 1", "-a e"), each = 2), length.out = N)
    tree <- rep("", N)
    if (!is.null(itree)) tree[] <- paste("-u ", itree)

    cmd <- paste(execname, "-i", seqfile, fmt, boot, mdl, tstv, inv, alpha, tree, "--append ")
    outfile <- paste(seqfile, "_phyml_stats.txt", sep = "")
    if (!append) {
        unlink(outfile)
        unlink(paste(seqfile, "_phyml_tree.txt", sep = ""))
    }
    imod <- 1:N
    if (!is.null(exclude)) imod <- imod[!.phymltest.model %in% exclude]

    for (i in imod) system(cmd[i])

    l <- readLines(outfile)
    l <- grep("Log-likelihood:", l, value = TRUE)
    ## in case there were already some results in the output file:
    if (dd <- length(l) - length(imod)) l <- l[-(1:dd)]
    loglik <- as.numeric(sub(". Log-likelihood:", "", l))
    names(loglik) <- .phymltest.model[imod]
    class(loglik) <- "phymltest"
    loglik
}

print.phymltest <- function(x, ...)
{
    nfp <- .phymltest.nfp[.phymltest.model %in% names(x)]
    X <- cbind(nfp, x, 2 * (nfp - x))
    rownames(X) <- names(x)
    colnames(X) <- c("nb.free.para", "loglik", "AIC")
    print(X)
}

summary.phymltest <- function(object, ...)
{
    nfp <- .phymltest.nfp[.phymltest.model %in% names(object)]
    N <- length(object)
    model1 <- model2 <- character(0)
    chi2 <- df <- P.val <- numeric(0)
    for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
            if (nfp[i] >= nfp[j]) next
            m1 <- unlist(strsplit(names(object)[i], "\\+"))
            m2 <- unlist(strsplit(names(object)[j], "\\+"))
            if (m1[1] == "K80" && m2[1] == "F81") next
            ## a verifier que ds les 2 lignes suivantes les conversions
            ## se font bien correctement!!!!
            if (length(grep("\\+I", names(object)[i])) > 0 && length(grep("\\+I", names(object)[j])) == 0) next
            if (length(grep("\\+G", names(object)[i])) > 0 && length(grep("\\+G", names(object)[j])) == 0) next
            ## Now we should be sure that m1 is nested in m2.
            chi2 <- c(chi2, 2 * (object[j] - object[i]))
            df <- c(df, nfp[j] - nfp[i])
            P.val <- c(P.val, 1 - pchisq(2 * (object[j] - object[i]), nfp[j] - nfp[i]))
            model1 <- c(model1, names(object)[i])
            model2 <- c(model2, names(object)[j])
        }
    }
    data.frame(model1, model2, chi2, df, P.val = round(P.val, 4))
}

plot.phymltest <- function(x, main = NULL, col = "blue", ...)
{
    nfp <- .phymltest.nfp[.phymltest.model %in% names(x)]
    N <- length(x)
    aic <- 2 * (nfp - x)
    if (is.null(main))
      main <- paste("Akaike information criterion for",
                    deparse(substitute(x)))
    plot(rep(1, N), aic, bty = "n", xaxt = "n", yaxt = "n",
         type = "n", xlab = "", ylab = "", main = main, ...)
    axis(side = 2, pos = 0.85, las = 2)
    abline(v = 0.85)
    y.lab <- seq(min(aic), max(aic), length = N)
    segments(0.85, sort(aic), 1.1, y.lab, col = col)
    text(1.1, y.lab,
         parse(text = sub("\\+G", "\\+Gamma", names(sort(aic)))),
         adj = 0)
}
