'print.parafit' <-
    function(x, ...)
{
	cat("\nTest of host-parasite coevolution",'\n','\n')

	cat("Global test:  ParaFitGlobal =",x$ParaFitGlobal,", p-value =", x$p.global, "(", x$nperm,"permutations)",'\n','\n')

	n.links <- nrow(x$link.table)
	cat("There are",n.links,"host-parasite links in matrix HP",'\n','\n')
	cat("Test of individual host-parasite links", "(", x$nperm, "permutations)",'\n','\n')
	print(x$link.table)

	cat('\n',"Number of parasites per host",'\n')
	print(x$para.per.host)
	cat('\n',"Number of hosts per parasite",'\n')
	print(x$host.per.para)
    invisible(x)
}
