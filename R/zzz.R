# $Id: zzz.R,v 1.8 2003/05/30 15:52:19 warnesgr Exp $

.First.lib <- function(libname, pkgname)
{
  if( is.null(formals(library)$pos) )
    library <- library.pos
  
  library(combinat,pos=3)
  library(MASS, pos=3)
  library(gregmisc,pos=3)
}
