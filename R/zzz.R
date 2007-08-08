## $Id: zzz.R 1297 2007-08-08 14:51:48Z warnes $

# Obsoleted by Proper use the DESCRIPTION and NAMESPACE files
.onLoad <- .First.lib <- function(libname, pkgname)
{
	cat("\n")
        cat("NOTE:\n")
	cat("\n")
	cat("  The R-Genetics project has developed an set of enhanced genetics\n")
	cat("  packages that will shortly replace 'genetics'. Please visit the \n")
        cat("  project homepage at http://rgenetics.org for more information.\n")
	cat("\n")
	
}
