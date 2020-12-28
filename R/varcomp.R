## varcomp.R (2004-10-29)

##   Variance Component of Mixed-Effect Linear Model

## Copyright 2004 Julien Dutheil

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

varcomp <- function(x, scale = FALSE, cum = FALSE)
{
  if (!("lme" %in% class(x))) stop("Object \"x\" is not of class \"lme\"")
  res <- seq(along = x$modelStruct$reStruct)
  var <- vector(length = length(res) + 1)
  for(i in res) {
    var[length(var) - i] <- attr(summary(x$modelStruct$reStruct[[i]]),"stdDev")[1]*x$sigma
  }
  var[length(var)] <- x$sigma
  var <- var^2
  if(scale) var <- var/sum(var)
  if(cum) var <- cumsum(var)
  names(var) <- c(rev(names(x$modelStruct$reStruct)), "Within")
  class(var) <- "varcomp"
  return(var)
}

plot.varcomp <- function(x, xlab = "Levels", ylab = "Variance", type = "b", ...) {
  if (!("varcomp" %in% class(x))) stop("Object \"x\" is not of class \"varcomp\"")
  return(xyplot(x ~ ordered(names(x), levels=rev(names(x))), xlab=xlab, ylab=ylab, type=type, ...))
}

# For debuging:
#data(carnivora)
#m <- lme(log10(SW) ~ 1, random = ~ 1|Order/SuperFamily/Family/Genus, data=carnivora)
#v <- varcomp(m,T,T)
#plot(v)

