# $Id: summary.genotype.R,v 1.6 2002/11/08 21:07:28 warnesgr Exp $
#
# $Log: summary.genotype.R,v $
# Revision 1.6  2002/11/08 21:07:28  warnesgr
#
# - DESCRIPTION: Updated version number and date
# - TODO: Updated todo list.
#
# Revision 1.5  2002/10/28 18:19:19  warnesgr
#
# - Fixed for syntax error in summary.genotype()
#
# Revision 1.4  2002/06/18 19:38:41  warnesgr
#
# Changes to fix problems reported by R CMD check.
#
# Revision 1.3  2001/06/13 20:00:22  warnes
# - Added "allele.names" attribute.
#
# Revision 1.2  2001/05/29 10:27:08  leisch
# add NA line only if NA's are present
#
# Revision 1.1  2001/05/25 23:37:06  warnes
#
# Moved summary.genotype and print.summary.genotype into a separate file
# (R/summary.genotype.R).  Added man/summary.genotype.Rd containing
# documentation for summary.genotype.
#
#

###
### Provide the frequency and proportions of alleles and genotypes
###
### Note: specifying any value for the parameter maxsum will cause fallback to
###       summary.factor.  This is so that summary.dataframe will give
###       reasonable output when it contains a genotype column.
###       (Hopefully we can figure out something better to do in this
###       case.)

summary.genotype  <-  function(object,...,maxobjectsum)
  {
    # if we are called from within summary.data.frame, fall back to
    # summary.factor so that we don't mess up the display
    if(!missing(maxobjectsum))
      return(NeobjecttMethod("summary",object,...,maxobjectsum))
    
    retval  <-  list()
#    retval$genotype  <- object
    retval$allele.names  <- allele.names(object)
    
    retval$locus  <- attr(object,"locus")
    class(retval)  <- "summary.genotype"
    af  <- table(allele(object))
    retval$allele.freq    <- cbind("Count"=af,"Proportion"=prop.table(af))

    gf  <- table( object )
    retval$genotype.freq  <- cbind("Count"=gf,"Proportion"=prop.table(gf))
    
    if(any(is.na(object))){
        retval$allele.freq    <- rbind(retval$allele.freq,
                                       "NA"=c(sum(is.na(allele(object))),NA))
        retval$genotype.freq  <- rbind(retval$genotype.freq,
                                       "NA"=c(sum(is.na(object)),NA))
    }
    
    return(retval)
  }

print.summary.genotype  <-  function(x,...,round=2)
  { 
    if(!is.null(x$locus))
      {
        cat("\n")
        print( x$locus )
      }
     
    cat("\n")
    cat("Allele Frequency:\n")
    print(round(x$allele.freq,digits=round),...)

    cat("\n")
    cat("Genotype Frequency:\n")
    print(round(x$genotype.freq,digits=round),...)
    
    cat("\n")

    invisible(x)
  }

