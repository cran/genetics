# $Id: zzz.R,v 1.3 2002/11/12 17:57:54 warnesgr Exp $
#
# $Log: zzz.R,v $
# Revision 1.3  2002/11/12 17:57:54  warnesgr
# - diseq.ci now requires the combinat library instead of boot.
#
# Revision 1.2  2002/10/11 18:01:00  warnesgr
#
# Updated version number.
#

.First.lib <- function(libname, pkgname)
  {

    if(!require(combinat))
      warning("Unable to load 'combinat' library.  Function `diseq.ci' will fail.")
  }
