## reconstruct.R (2016-07-15)

##   Ancestral Character Estimation

## Copyright 2014-2016 Manuela Royer-Carenzi, Didier Gilles

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

#renvoie la racine d'arbre
racine <- function(arbre) {
    Ntip(arbre) + 1
}



# renvoie une liste dont la premiere composante est l'arbre renumerote
# de telle sorte que l'index d'un enfant est superieur a celui de son pere,
# la seconde compopsante est la fonction de l'index initial vers le second,
# et la troisieme son inverse
# (attention probleme pour l'image de 0 mise a l'image du max)
#
renumeroteArbre <- function(arbre) {
  m <- Ntip(arbre) + Nnode(arbre)
  v<-numeric(m)
  t<-numeric(m)
  stack<-numeric(m)
  istack<-1
  stack[istack]<-racine(arbre)
  codeI<-1
  codeL<-Nnode(arbre)+1
  while(istack>0){
    cour<-stack[istack]
    istack<-istack-1
    l <- which(arbre$edge[, 1] == cour)
    if(length(l)>0){
      v[cour] <- codeI
      t[codeI] <- cour
      codeI <- codeI+1
      for(i in 1:length(l)) {
	istack<-istack+1
	stack[istack] <- arbre$edge[l[i], 2]
      }
    } else {
      v[cour] <- codeL
      t[codeL] <- cour
      codeL <- codeL+1
    }
  }
  arbrebis<-arbre
#renumeroter les noms
  for(i in 1:Nedge(arbre)) {
    arbrebis$edge[i,1] <- v[arbre$edge[i,1]]
    arbrebis$edge[i,2] <- v[arbre$edge[i,2]]
  }
  l <- list(arbre = arbrebis, cod = v, dec = t)
  l
}


#calcule la matrice C selon le modele BM ou ABM
#
calculeC_ABM <- function(arbre) {
  m <- max(arbre[["edge"]])
  C <- matrix(0,nrow=m,ncol=m)
  for(i in 1:(m)) {
    l <- which(arbre$edge[, 2] == i)
    if(length(l)>0){
      for(j in 1:(m)) {
	C[j,i] <- C[j, arbre$edge[l[1], 1]]
      }
    }
    C[i,i]<-1;
  }
  t(C)
}

#calcule la matrice C selon le modele OU ou OU*
#
calculeC_OU <- function(arbre, a) {
  m <- max(arbre[["edge"]])
  C <- matrix(0,nrow=m,ncol=m)
  for(i in 1:(m)) {
    l <- which(arbre$edge[, 2] == i)
    if(length(l)>0){
      for(j in 1:(m)) {
	  C[j,i] <- C[j, arbre$edge[l[1], 1]]*exp(-a*arbre$edge.length[l[1]])
      }
    }
    C[i,i]<-1;
  }
  t(C)
}

#calcule la matrice C selon le modele type qui vaut ABM ou OU
calculeC <- function(type, arbre, a) {
  switch(type, ABM = calculeC_ABM(arbre), OU = calculeC_OU(arbre, a))
}





#########################
#calcul Variance
#########################

getSumSquare <- function(value, arbre) {
 sum <- 0.
 for(eu in 1:Nedge(arbre))
   sum <- sum + (value[arbre$edge[eu,2]]-value[arbre$edge[eu,1]])^2/arbre$edge.length[eu]
 sum
}


getMLHessian <- function(value, arbre) {
   sumSqu <- getSumSquare(value, arbre)
   nI <- Nnode(arbre)
   nT <- length(arbre$tip.label)
   nE <- nI+nT-1
   sizeH<-nI+1
   hessian <- matrix(0., nrow=sizeH, ncol=sizeH)
   var <- sumSqu/nE
   sd <- sqrt(var)
   hessian[1,1] <- -nE/(2*var^2)+sumSqu/var^3
   for(i in 1:nI) {
     	child <- which(arbre$edge[, 1] == nT+i)
  	if(length(child)>0) {
      		for(j in 1:length(child)) {
     		hessian[1,i+1] <- hessian[1,i+1]-(value[arbre$edge[child[j],2]]-value[nT+i])/arbre$edge.length[child[j]]
      		hessian[i+1,i+1] <- hessian[i+1,i+1]+1./arbre$edge.length[child[j]]
       		if(arbre$edge[child[j],2]>nT)
		  hessian[i+1,arbre$edge[child[j],2]-nT+1] <- -1./(var*arbre$edge.length[child[j]])
      		}
     	}
      	anc <- which(arbre$edge[, 2] == nT+i)
     	if(length(anc)>0) {
     	 	for(j in 1:length(anc)) {
       		hessian[1,i+1] <- hessian[1,i+1]+(value[nT+i]-value[arbre$edge[anc[j],1]])/arbre$edge.length[anc[j]]
     		hessian[i+1,i+1] <- hessian[i+1,i+1]+1./arbre$edge.length[anc[j]]
       		hessian[i+1,arbre$edge[anc[j],1]-nT+1] <- -1./(var*arbre$edge.length[anc[j]])
      		}
     	}
   hessian[1,i+1] <- -hessian[1,i+1]/sd^2
   hessian[i+1,1] <- hessian[1,i+1]
   hessian[i+1,i+1] <- hessian[i+1,i+1]/var
   }
   hessian
}

getREMLHessian <- function(value, arbre, sigma2) {
   nI <- Nnode(arbre)
   nT <- length(arbre$tip.label)
   sizeH<-nI
   hessian <- matrix(0., nrow=sizeH, ncol=sizeH)
   for(i in 1:nI) {
     	child <- which(arbre$edge[, 1] == nT+i)
  	if(length(child)>0) {
      		for(j in 1:length(child)) {
      		hessian[i,i] <- hessian[i,i]+1./arbre$edge.length[child[j]]
       		if(arbre$edge[child[j],2]>nT)
		  hessian[i,arbre$edge[child[j],2]-nT] <- -1./(sigma2*arbre$edge.length[child[j]])
      		}
     	}
      	anc <- which(arbre$edge[, 2] == nT+i)
     	if(length(anc)>0) {
     	 	for(j in 1:length(anc)) {
      		hessian[i,i] <- hessian[i,i]+1./arbre$edge.length[anc[j]]
       		hessian[i,arbre$edge[anc[j],1]-nT] <- -1./(sigma2*arbre$edge.length[anc[j]])
      		}
     	}
      hessian[i,i] <- hessian[i,i]/sigma2
   }
   hessian
}



glsBM <- function (phy, x, CI=TRUE) 
{	
	obj <- list()
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nbTotN <- nb.tip+nb.node
	sigmaMF <- 1	
	TsTemps <- dist.nodes(phy)
	TempsRacine <- TsTemps[(nb.tip+1),]
	IndicesMRCA <- mrca(phy, full=T)
	M <- matrix(NA, ncol=nbTotN, nrow=nbTotN)
	for (i in 1:nbTotN)
	{
	for (j in 1:nbTotN)
	{
		M[i,j] <- sigmaMF^2 * TempsRacine[IndicesMRCA[i,j]]
	}
	}
# M = SigmaZ
	varAL <- M[-(1:nb.tip), 1:nb.tip]
	varAA <- M[-(1:nb.tip), -(1:nb.tip)]
        varLL <- M[(1:nb.tip), 1:nb.tip]
        invVarLL <- solve(varLL)
	UL <- rep(1, length=nb.tip)
	UA <- rep(1, length=nb.node)
	TL <- TempsRacine[1:nb.tip] 
	TA <- TempsRacine[(nb.tip+1):(nb.tip+nb.node)]
#
	IVL_Z <- invVarLL %*% x
	IVL_T <- invVarLL %*% TL
	IVL_U <- invVarLL %*% UL
	U_IVL_U <- t(UL) %*% IVL_U
	U_IVL_Z <- t(UL) %*% IVL_Z
	DeltaU <- UA - varAL %*% IVL_U
#
	Racine_chap <- solve(U_IVL_U) %*% U_IVL_Z 
	Racine_chap <- as.numeric(Racine_chap)
	Anc_chap <- Racine_chap * DeltaU + varAL %*% IVL_Z
	Anc_chap <- as.vector(Anc_chap)
	obj$ace <- Anc_chap
	names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
#
 	if (CI) {
		Vec <- x - Racine_chap
		Num <- t(Vec) %*% invVarLL %*% Vec
		Num <- as.numeric(Num)
		sigma2_chap <- Num / (nb.tip-1)
		obj$sigma2 <- sigma2_chap
                se <- sqrt((varAA - varAL %*% invVarLL %*% t(varAL))[cbind(1:nb.node,
                  1:nb.node)])
		se <- se * sqrt(sigma2_chap)
                tmp <- se * qt(0.025, df=nb.tip-1)
                 obj$CI95 <- cbind(lower=obj$ace + tmp, upper=obj$ace - tmp)
           }
	obj
}

glsABM <- function (phy, x, CI=TRUE)
{
	obj <- list()
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nbTotN <- nb.tip+nb.node
	sigmaMF <- 1
	TsTemps <- dist.nodes(phy)
	TempsRacine <- TsTemps[(nb.tip+1),]
	IndicesMRCA <- mrca(phy, full=T)
	M <- matrix(NA, ncol=nbTotN, nrow=nbTotN)
	for (i in 1:nbTotN)
	{
	for (j in 1:nbTotN)
	{
		M[i,j] <- sigmaMF^2 * TempsRacine[IndicesMRCA[i,j]]
	}
	}
# M = SigmaZ
	varAL <- M[-(1:nb.tip), 1:nb.tip]
	varAA <- M[-(1:nb.tip), -(1:nb.tip)]
        varLL <- M[(1:nb.tip), 1:nb.tip]
        invVarLL <- solve(varLL)
	UL <- rep(1, length=nb.tip)
	UA <- rep(1, length=nb.node)
	TL <- TempsRacine[1:nb.tip]
	TA <- TempsRacine[(nb.tip+1):(nb.tip+nb.node)]
#
	IVL_Z <- invVarLL %*% x
	IVL_T <- invVarLL %*% TL
	IVL_U <- invVarLL %*% UL
	U_IVL_U <- t(UL) %*% IVL_U
	T_IVL_T <- t(TL) %*% IVL_T
	U_IVL_T <- t(UL) %*% IVL_T
	U_IVL_Z <- t(UL) %*% IVL_Z
	T_IVL_Z <- t(TL) %*% IVL_Z
	DeltaT <- TA - varAL %*% IVL_T
	DeltaU <- UA - varAL %*% IVL_U
#
	Den <- U_IVL_U * T_IVL_T - U_IVL_T^2
	Den <- as.numeric(Den)
	Mu_chap <- (U_IVL_U * T_IVL_Z - U_IVL_T * U_IVL_Z) / Den
	Mu_chap <- as.numeric(Mu_chap)
	Racine_chap <- (T_IVL_T * U_IVL_Z - U_IVL_T * T_IVL_Z) / Den
	Racine_chap <- as.numeric(Racine_chap)
	Anc_chap <- Mu_chap * DeltaT + Racine_chap * DeltaU + varAL %*% IVL_Z
	Anc_chap <- as.vector(Anc_chap)
	obj$ace <- Anc_chap
	names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
	obj$mu <- Mu_chap
#
 	if (CI) {
		Vec <- x - Racine_chap - Mu_chap * TL
		Num <- t(Vec) %*% invVarLL %*% Vec
		Num <- as.numeric(Num)
		sigma2_chap <- Num / (nb.tip-2)
		obj$sigma2 <- sigma2_chap
                se <- sqrt((varAA - varAL %*% invVarLL %*% t(varAL))[cbind(1:nb.node, 
                  1:nb.node)])
		se <- se * sqrt(sigma2_chap)
                tmp <- se * qt(0.025, df=nb.tip-2)
                obj$CI95 <- cbind(lower=obj$ace + tmp, upper=obj$ace - tmp)
            }
	obj
}

# theta = z0
glsOUv1 <- function (phy, x, alpha, CI=TRUE)
{
	obj <- list()
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nbTotN <- nb.tip+nb.node
	sigmaMF <- 1
	alphaM <- alpha
	nbTotN <- nb.tip+nb.node
	TsTemps <- dist.nodes(phy)
	TempsRacine <- TsTemps[(nb.tip+1),]
	IndicesMRCA <- mrca(phy, full=T)
	M <- matrix(NA, ncol=nbTotN, nrow=nbTotN)
	for (i in 1:nbTotN)
		{
	for (j in 1:nbTotN)
		{
		Tempsm <- TempsRacine[IndicesMRCA[i,j]]
		Tempsi <- TempsRacine[i]
		Tempsj <- TempsRacine[j]
		M[i,j] <- sigmaMF^2 * exp(-alphaM * (Tempsi+Tempsj-2*Tempsm)) * (1-exp(-2*alphaM * Tempsm)) / (2 * alphaM)
		}
		}
# M = SigmaZ
	varAL <- M[-(1:nb.tip), 1:nb.tip]
	varAA <- M[-(1:nb.tip), -(1:nb.tip)]
        varLL <- M[(1:nb.tip), 1:nb.tip]
        invVarLL <- solve(varLL)
	UL <- rep(1, length=nb.tip)
	UA <- rep(1, length=nb.node)
	TL <- TempsRacine[1:nb.tip]
	TA <- TempsRacine[(nb.tip+1):(nb.tip+nb.node)]
#
	IVL_Z <- invVarLL %*% x
	IVL_T <- invVarLL %*% TL
	IVL_U <- invVarLL %*% UL
	U_IVL_U <- t(UL) %*% IVL_U
	U_IVL_Z <- t(UL) %*% IVL_Z
	DeltaU <- UA - varAL %*% IVL_U
#
	Racine_chap <- solve(U_IVL_U) %*% U_IVL_Z
	Racine_chap <- as.numeric(Racine_chap)
	Anc_chap <- Racine_chap * DeltaU + varAL %*% IVL_Z
	Anc_chap <- as.vector(Anc_chap)
	obj$ace <- Anc_chap
	names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
#
# vraisemblance
#
	mL <- Racine_chap
	Num <- t(x-mL) %*% invVarLL %*% (x-mL)
	Num <- as.numeric(Num)
	sigma2_chap <- Num / (nb.tip-1)
	obj$sigma <- sqrt(sigma2_chap)
	VL <- sigma2_chap * varLL
	invVL <- invVarLL / sigma2_chap
	LVrais <- - t(x-mL) %*% invVL %*% (x-mL) /2 - nb.tip * log(2*pi)/2 - log(det(VL))/2
	obj$LLik <- as.numeric(LVrais)
#
 	if (CI) {
                se <- sqrt((varAA - varAL %*% invVarLL %*% t(varAL))[cbind(1:nb.node,
                  1:nb.node)])
		se <- se * sqrt(sigma2_chap)
                tmp <- se * qt(0.025, df=nb.tip-1)
                 obj$CI95 <- cbind(lower=obj$ace + tmp, upper=obj$ace - tmp)
            }
	obj
}


# theta pas egal a z0
glsOUv2 <- function (phy, x, alpha, CI=TRUE)
{
	obj <- list()
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nbTotN <- nb.tip+nb.node
	sigmaMF <- 1
	nbTotN <- nb.tip+nb.node
	TsTemps <- dist.nodes(phy)
	TempsRacine <- TsTemps[(nb.tip+1),]
	IndicesMRCA <- mrca(phy, full=T)
	M <- matrix(NA, ncol=nbTotN, nrow=nbTotN)
	for (i in 1:nbTotN)
		{
	for (j in 1:nbTotN)
		{
		Tempsm <- TempsRacine[IndicesMRCA[i,j]]
		Tempsi <- TempsRacine[i]
		Tempsj <- TempsRacine[j]
		M[i,j] <- sigmaMF^2 * exp(-alpha * (Tempsi+Tempsj-2*Tempsm)) * (1-exp(-2*alpha * Tempsm)) / (2 * alpha)
		}
		}
# M = SigmaZ
	varAL <- M[-(1:nb.tip), 1:nb.tip]
	varAA <- M[-(1:nb.tip), -(1:nb.tip)]
        varLL <- M[(1:nb.tip), 1:nb.tip]
        invVarLL <- solve(varLL)
	vecW <- exp(-alpha * TempsRacine)
	UL <- vecW[1:nb.tip]
	UA <- vecW[(nb.tip+1):(nb.tip+nb.node)]
	TL <- 1-UL
	TA <- 1-UA
#
#
	IVL_Z <- invVarLL %*% x
	IVL_T <- invVarLL %*% TL
	IVL_U <- invVarLL %*% UL
	U_IVL_U <- t(UL) %*% IVL_U
	T_IVL_T <- t(TL) %*% IVL_T
	U_IVL_T <- t(UL) %*% IVL_T
	U_IVL_Z <- t(UL) %*% IVL_Z
	T_IVL_Z <- t(TL) %*% IVL_Z
	DeltaT <- TA - varAL %*% IVL_T
	DeltaU <- UA - varAL %*% IVL_U
#
	Den <- U_IVL_U * T_IVL_T - U_IVL_T^2
	Den <- as.numeric(Den)
	Theta_chap <- (U_IVL_U * T_IVL_Z - U_IVL_T * U_IVL_Z) / Den
	Theta_chap <- as.numeric(Theta_chap)
	Racine_chap <- (T_IVL_T * U_IVL_Z - U_IVL_T * T_IVL_Z) / Den
	Racine_chap <- as.numeric(Racine_chap)
	Anc_chap <- Theta_chap * DeltaT + Racine_chap * DeltaU + varAL %*% IVL_Z
	Anc_chap <- as.vector(Anc_chap)
	obj$ace <- Anc_chap
	names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
	obj$theta <- Theta_chap
#
# vraisemblance
#
	mL <- (Racine_chap * UL + Theta_chap * TL)
	Num <- t(x-mL) %*% invVarLL %*% (x-mL)
	Num <- as.numeric(Num)
	sigma2_chap <- Num / (nb.tip-2)
	obj$sigma <- sqrt(sigma2_chap)
	VL <- sigma2_chap * varLL
	invVL <- invVarLL / sigma2_chap
	LVrais <- - t(x-mL) %*% invVL %*% (x-mL) /2 - nb.tip * log(2*pi)/2 - log(det(VL))/2
	obj$LLik <- as.numeric(LVrais)
#
 	if (CI) {
                se <- sqrt((varAA - varAL %*% invVarLL %*% t(varAL))[cbind(1:nb.node,
                  1:nb.node)])
		se <- se * sqrt(sigma2_chap)
                tmp <- se * qt(0.025, df=nb.tip-2)
                obj$CI95 <- cbind(lower=obj$ace + tmp, upper=obj$ace - tmp)
            }
	obj
}


reconstruct <- function (x, phyInit, method = "ML", alpha = NULL, CI = TRUE) {
 if(!is.null(alpha)) {
  if(alpha<=0)
   stop("alpha has to be positive.")
 }
 if (!inherits(phyInit, "phylo"))
  stop("object \"phy\" is not of class \"phylo\"")
 if (is.null(phyInit$edge.length))
  stop("tree has no branch lengths")
 nN <- phyInit$Nnode
 nT <- length(x)
 switch(method,
  ML = {
   Intern <- glsBM(phy=phyInit, x=x, CI=F)$ace
   Value <- c(x, Intern)
   Hessian <- getMLHessian(Value, phyInit)
   se <- sqrt(diag(solve(Hessian)))
   se <- se[-1]
   tmp <- se*qt(0.025, nN)
   CI95 <- cbind(lower=Intern+tmp, upper=Intern-tmp)
  },
  REML={
   minusLogLik <- function(sig2) {
    if (sig2 < 0) return(1e+100)
    V <- sig2 * vcv(phyInit)
    distval <- stats::mahalanobis(x, center = mu, cov = V)
    logdet <- sum(log(eigen(V, symmetric = TRUE, only.values = TRUE)$values))
    (nT * log(2 * pi) + logdet + distval)/2
   }
   Intern <- glsBM(phy=phyInit, x=x, CI=F)$ace
   Value <- c(x, Intern)
   GM <- Intern[1]
   mu <- rep(GM, nT)
   out <- nlm(minusLogLik, 1, hessian = FALSE)
   sigma2 <- out$estimate
   Hessian <- getREMLHessian(Value, phyInit, sigma2)
   se <- sqrt(diag(solve(Hessian)))
   tmp <- se*qt(0.025, nN)
   CI95 <- cbind(lower=Intern+tmp, upper=Intern-tmp)
  },
  GLS = {
	result <- glsBM(phy=phyInit, x=x, CI=T)
    Intern <- result$ace
    CI95 <- result$CI95
  },
  GLS_ABM = {
	result <- glsABM(phy=phyInit, x=x, CI=T)
    Intern <- result$ace
    CI95 <- result$CI95
  },
  GLS_OUS = {
	if(is.null(alpha)) {
		funOpt1 <- function(alpha)
		{
			-glsOUv1(phy=phyInit, x=x, alpha, CI=F)$LLik
		}
		calOp <- optim(par=0.25, fn=funOpt1, method="Brent", lower=0.001, upper=1)
		if (calOp$convergence == 0)
		{
			alpha <- calOp$par
		} else {
			stop("Estimation error for alpha")
		}
	}
	result <- glsOUv1(phy=phyInit, x=x, alpha=alpha, CI=T)
    Intern <- result$ace
    CI95 <- result$CI95
  },
  GLS_OU = {
	if(is.null(alpha)) {
		funOpt2 <- function(alpha)
		{
			-glsOUv2(phy=phyInit, x=x, alpha, CI=F)$LLik
		}
		calOp <- optim(par=0.25, fn=funOpt2, method="Brent", lower=0.001, upper=1)
		if (calOp$convergence == 0)
		{
			alpha <- calOp$par
		} else {
			stop("Estimation error for alpha")
		}
	}
	result <- glsOUv2(phy=phyInit, x=x, alpha=alpha, CI=T)
    Intern <- result$ace
    CI95 <- result$CI95
  }
 )
if (CI==TRUE)
 list(ace=Intern, CI95=CI95)
else
  list(ace=Intern)
}
