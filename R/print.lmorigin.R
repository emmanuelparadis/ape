'print.lmorigin' <-
    function(x, ...)
{
	if(x$origin) {
		cat("\nRegression through the origin",'\n')
		} else {
		cat("\nMultiple regression with estimation of intercept",'\n')
		}
	cat("\nCall:\n")
    cat(deparse(x$call),'\n')
	if(x$origin) { names <- x$var.names[-1] }
		else { names <- c("(Intercept)",x$var.names[-1]) }

	cat("\nCoefficients and parametric test results \n",'\n')
	res <- as.data.frame(cbind(summary(x$reg)$coefficients[,1], summary(x$reg)$coefficients[,2], summary(x$reg)$coefficients[,3], summary(x$reg)$coefficients[,4]))
	rownames(res) <- names
	colnames(res) <- c("Coefficient","Std_error","t-value","Pr(>|t|)")
	printCoefmat(res, P.values=TRUE, signif.stars=TRUE)
	
	if(x$nperm > 0) {
		cat("\nTwo-tailed tests of regression coefficients\n",'\n')
		res2 <- as.data.frame(cbind(summary(x$reg)$coefficients[,1], x$p.param.t.2tail, x$p.perm.t.2tail))
		rownames(res2) <- names
		colnames(res2) <- c("Coefficient","p-param","p-perm")
		nc <- 3
		printCoefmat(res2, P.values=TRUE, signif.stars=TRUE, has.Pvalue = 3 && substr(colnames(res2)[3],1,6) == "p-perm")

		cat("\nOne-tailed tests of regression coefficients:",'\n')
		cat("test in the direction of the sign of the coefficient\n",'\n')
		res1 <- as.data.frame(cbind(summary(x$reg)$coefficients[,1], x$p.param.t.1tail, x$p.perm.t.1tail))
		rownames(res1) <- names
		colnames(res1) <- c("Coefficient","p-param","p-perm")
		nc <- 3
		printCoefmat(res1, P.values=TRUE, signif.stars=TRUE, has.Pvalue = 3 && substr(colnames(res1)[3],1,6) == "p-perm")

		}
	cat("\nResidual standard error:", summary(x$reg)$sigma, "on", summary(x$reg)$df[2],"degrees of freedom",'\n')
	cat("Multiple R-square:", summary(x$reg)$r.squared,"  Adjusted R-square:", summary(x$reg)$adj.r.squared,'\n')

	F   <- summary(x$reg)$fstatistic[[1]]
	df1 <- summary(x$reg)$fstatistic[[2]]
	df2 <- summary(x$reg)$fstatistic[[3]]
	p.param.F <- pf(F, df1, df2, lower.tail=FALSE)
	cat("\nF-statistic:", F, "on", df1, "and", df2, "DF:\n")
	cat("   parametric p-value   :", p.param.F,'\n')
	if(x$nperm > 0) {
		cat("   permutational p-value:", x$p.perm.F,'\n')
		if(x$method == "raw") {
			cat("after",x$nperm,"permutations of",x$method,"data",'\n','\n')
			} else {
			cat("after",x$nperm,"permutations of",x$method,"of full model",'\n','\n')
			}			
		}
    invisible(x) 
}
