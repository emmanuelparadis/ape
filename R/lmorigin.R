'lmorigin' <- 
function(formula, data=NULL, origin=TRUE, nperm=999, method=NULL, silent=FALSE)
#
# This program computes a multiple linear regression and performs tests 
# of significance of the equation parameters using permutations. 
# 
# origin=TRUE: the regression line can be forced through the origin. Testing 
# the significance in that case requires a special permutation procedure.
#
# Permutation methods: raw data or residuals of full model
#    Default method in regression through the origin: raw data
#    Default method in ordinary multiple regression: residuals of full model
#    - In ordinary multiple regression when m = 1: raw data
#
#         Pierre Legendre, March 2009
{
if(!is.null(method)) method <- match.arg(method, c("raw", "residuals"))
if(is.null(method) & origin==TRUE) method <- "raw"
if(is.null(method) & origin==FALSE) method <- "residuals"
if(nperm < 0) stop("Incorrect value for 'nperm'")
	
## From the formula, find the variables and the number of observations 'n'
toto <- lm(formula, data)
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "offset"), names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
var.names = colnames(mf)   # Noms des variables
y <- as.matrix(mf[,1])
colnames(y) <- var.names[1]
X <- as.matrix(mf[,-1])
n <- nrow(mf)
m <- ncol(X)

a <- system.time({
mm<- m           # No. regression coefficients, possibly including the intercept
if(m == 1) method <- "raw"
if(nrow(X) != n) stop("Unequal number of rows in y and X")

if(origin) {
	if(!silent) cat("Regression through the origin",'\n')
	reg <- lm(y ~ 0 + X)
	} else {
	if(!silent) cat("Multiple regression with estimation of intercept",'\n')
	reg <- lm(y ~ X)	
	mm <- mm+1
	}

if(!silent) {
	if(nperm > 0) {
		if(method == "raw") {
			cat("Permutation method =",method,"data",'\n')
			} else {
			cat("Permutation method =",method,"of full model",'\n')
			}
		}
	}
t.vec      <- summary(reg)$coefficients[,3]
p.param.t  <- summary(reg)$coefficients[,4]
df1        <- summary(reg)$fstatistic[[2]]
df2        <- summary(reg)$fstatistic[[3]]
F          <- summary(reg)$fstatistic[[1]]
y.res      <- summary(reg)$residuals
# b.vec      <- summary(reg)$coefficients[,1]
# r.sq       <- summary(reg)$r.squared
# adj.r.sq   <- summary(reg)$adj.r.squared
# p.param.F  <- pf(F, df1, df2, lower.tail=FALSE)

if(df1 < m) stop("\nCollinearity among the X variables. Check using 'lm'")

# Permutation tests
if(nperm > 0) {

	nGT.F  <- 1
	nGT1.t <- rep(1,mm)
	nGT2.t <- rep(1,mm)
	sign.t <- sign(t.vec)

	for(i in 1:nperm)  # Permute raw data. Always use this method for F-test
		{
		if(origin) {   # Regression through the origin
			dia.bin <- diag((rbinom(n,1,0.5)*2)-1)
			y.perm <- dia.bin %*% sample(y)
			reg.perm <- lm(y.perm ~ 0 + X)
			} else {   # Multiple linear regression
			y.perm <- sample(y,n)
			reg.perm <- lm(y.perm ~ X)
			}
	
		# Permutation test of the F-statistic
		F.perm <- summary(reg.perm)$fstatistic[1]
		if(F.perm >= F) nGT.F <- nGT.F+1

		# Permutation tests of the t-statistics: permute raw data
		if(method == "raw") {
			t.perm <- summary(reg.perm)$coefficients[,3]
			if(nperm <= 5) cat(t.perm,'\n')
			for(j in 1:mm) {
				# One-tailed test in direction of sign
				if(t.perm[j]*sign.t[j] >= t.vec[j]*sign.t[j]) nGT1.t[j] <- nGT1.t[j]+1
				# Two-tailed test
				if( abs(t.perm[j]) >= abs(t.vec[j]) ) nGT2.t[j] <- nGT2.t[j]+1
				}
			}
		}
		
	if(method == "residuals") {
	# Permute residuals of full model
		for(i in 1:nperm) {
			if(origin) {   # Regression through the origin
				dia.bin <- diag((rbinom(n,1,0.5)*2)-1)
				y.perm <- dia.bin %*% sample(y.res)
				reg.perm <- lm(y.perm ~ 0 + X)
				} else {   # Multiple linear regression
				y.perm <- sample(y.res,n)
				reg.perm <- lm(y.perm ~ X)
				}
	
			# Permutation tests of the t-statistics: permute residuals
			t.perm <- summary(reg.perm)$coefficients[,3]
			if(nperm <= 5) cat(t.perm,'\n')
			for(j in 1:mm) {
				# One-tailed test in direction of sign
				if(t.perm[j]*sign.t[j] >= t.vec[j]*sign.t[j]) nGT1.t[j] <- nGT1.t[j]+1
				# Two-tailed test
				if( abs(t.perm[j]) >= abs(t.vec[j]) ) nGT2.t[j] <- nGT2.t[j]+1
				}
			}
		}
	# Compute the permutational probabilities
	p.perm.F <- nGT.F/(nperm+1)
	p.perm.t1 <- nGT1.t/(nperm+1)
	p.perm.t2 <- nGT2.t/(nperm+1)

	### Do not test intercept by permutation of residuals in multiple regression
	if(!origin & method=="residuals") {
		if(silent) {   # Note: silent==TRUE in simulation programs
			p.perm.t1[1] <- p.perm.t2[1] <- 1
			} else {
			p.perm.t1[1] <- p.perm.t2[1] <- NA
			}
		}
	}

})
a[3] <- sprintf("%2f",a[3])
if(!silent) cat("Computation time =",a[3]," sec",'\n')
#
if(nperm == 0) {

	out <- list(reg=reg, p.param.t.2tail=p.param.t, p.param.t.1tail=p.param.t/2, origin=origin, nperm=nperm, var.names=var.names, call=match.call())

	} else {

	out <- list(reg=reg, p.param.t.2tail=p.param.t, p.param.t.1tail=p.param.t/2, p.perm.t.2tail=p.perm.t2, p.perm.t.1tail=p.perm.t1, p.perm.F=p.perm.F, origin=origin, nperm=nperm, method=method, var.names=var.names, call=match.call())

	}
#
class(out) <- "lmorigin"
out
}
