## $Id: zzz.R 1330 2007-11-20 20:23:12Z warnes $

# Obsoleted by Proper use the DESCRIPTION and NAMESPACE files
.onLoad <- .First.lib <- function(libname, pkgname)
{
	cat("\n")
        cat("NOTE: THIS PACKAGE IS NOW OBSOLETE.\n")
	cat("\n")
	cat("  The R-Genetics project has developed an set of enhanced genetics\n")
	cat("  packages to replace 'genetics'. Please visit the project homepage\n")
        cat("  at http://rgenetics.org for informtion.\n")
	cat("\n")
	
}
