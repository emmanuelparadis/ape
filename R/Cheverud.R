## Cheverud.R (2004-10-29)

##    Cheverud's 1985 Autoregression Model

## Copyright 2004 Julien Dutheil

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

# This function is adapted from a MatLab code from
# Rholf, F. J. (2001) Comparative Methods for the Analysis of Continuous Variables: Geometric Interpretations.
# Evolution 55(11): 2143-2160
compar.cheverud <- function(y, W, tolerance=1e-6, gold.tol=1e-4)
{
  ## fix by Michael Phelan
  diag(W) <- 0 # ensure diagonal is zero
  ## end of fix
  y <- as.matrix(y)
  if(dim(y)[2] != 1) stop("Error: y must be a single column vector.")
  D <- solve(diag(apply(t(W),2,sum)))
  Wnorm <- D %*% W #Row normalize W matrix
  n <- dim(y)[1]
  m <- dim(y)[2]
  y <- y-matrix(rep(1, n)) %*% apply(y,2,mean) # Deviations from mean
  Wy <- Wnorm %*% y

  Wlam <- eigen(Wnorm)$values # eigenvalues of W

  # Find distinct eigenvalues
  sorted <- sort(Wlam)
  # Check real:
  for (ii in 1:n) {
    if(abs(Im(sorted[ii])) > 1e-12) {
      warning(paste("Complex eigenvalue coerced to real:", Im(sorted[ii])))
	  }
    sorted[ii] <- Re(sorted[ii]) # Remove imaginary part
  }
  sorted <- as.double(sorted)

	Distinct <- numeric(0)
  Distinct[1] <- -Inf
  Distinct[2] <- sorted[1]
  nDistinct <- 2
  for(ii in 2:n) {
    if(sorted[ii] - Distinct[nDistinct] > tolerance) {
      nDistinct <- nDistinct + 1
      Distinct[nDistinct] <- sorted[ii]
    }
  }

  # Search for minimum of LL

  likelihood <- function(rhohat) {
    DetProd <- 1
    for(j in 1:n) {
      prod <- 1 - rhohat * Wlam[j]
      DetProd <- DetProd * prod
    }
    absValDet <- abs(DetProd) #[abs to allow rho > 1]
    logDet <- log(absValDet)
    LL <- log(t(y) %*% y - 2 * rhohat * t(y) %*% Wy + rhohat * rhohat * t(Wy) %*% Wy) - logDet*2/n
    return(LL)
  }

  GoldenSearch <- function(ax, cx) {
    # Golden section search over the interval ax to cx
    # Return rhohat and likelihood value.
    r <- 0.61803399
    x0 <- ax
    x3 <- cx
    bx <- (ax + cx)/2
    if(abs(cx - bx) > abs(bx - ax)) {
      x1 <- bx
      x2 <- bx + (1-r)*(cx - bx)
    } else {
      x2 <- bx
      x1 <- bx - (1-r)*(bx - ax)
    }
    f1 <- likelihood(x1)
    f2 <- likelihood(x2)
    while(abs(x3 - x0) > gold.tol*(abs(x1) + abs(x2))) {
      if(f2 < f1) {
        x0 <- x1
        x1 <- x2
        x2 <- r * x1 + (1 - r) * x3
        f1 <- f2
        f2 <- likelihood(x2)
      } else {
        x3 <- x2
        x2 <- x1
        x1 <- r * x2 + (1 - r) * x0
        f2 <- f1
        f1 <- likelihood(x1)
      }
    }
    if(f1 < f2) {
      likelihood <- f1
      xmin <- x1
    } else {
      likelihood <- f2
      xmin <- x2
    }
    return(list(rho=xmin, LL=likelihood))
  }

  LL <- Inf
  for(ii in 2:(nDistinct -1)) {# Search between pairs of roots
    # [ constrain do not use positive roots < 1]
    ax <- 1/Distinct[ii]
    cx <- 1/Distinct[ii+1]
    GS <- GoldenSearch(ax, cx)
    if(GS$LL < LL) {
      LL <- GS$LL
      rho <- GS$rho
    }
  }
  # Compute residuals:
  res <- y - rho * Wy
  return(list(rhohat=rho, Wnorm=Wnorm, residuals=res))
}

#For debugging:
#W<- matrix(c(
#  0,1,1,2,0,0,0,0,
#  1,0,1,2,0,0,0,0,
#  1,1,0,2,0,0,0,0,
#  2,2,2,0,0,0,0,0,
#  0,0,0,0,0,1,1,2,
#  0,0,0,0,1,0,1,2,
#  0,0,0,0,1,1,0,2,
#  0,0,0,0,2,2,2,0
#),8)
#W <- 1/W
#W[W == Inf] <- 0
#y<-c(-0.12,0.36,-0.1,0.04,-0.15,0.29,-0.11,-0.06)
#compar.cheverud(y,W)
#
#y<-c(10,8,3,4)
#W <- matrix(c(1,1/6,1/6,1/6,1/6,1,1/2,1/2,1/6,1/2,1,1,1/6,1/2,1,1), 4)
#compar.cheverud(y,W)

