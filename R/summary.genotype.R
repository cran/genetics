# $Id: summary.genotype.R,v 1.8 2003/05/22 17:25:23 warnesgr Exp $

###
### Provide the frequency and proportions of alleles and genotypes
###

# used when summary.genotype is called from summary.data.frame:
shortsummary.genotype <- function(object, ..., maxsum)
  {
    tmp <- summary.factor(object, maxsum=maxsum)
    retval <- paste(format(tmp), " (", format(round(prop.table(tmp)*100)), "%)", sep='' )
    names(retval) <- names(tmp)
    #retval <- retval[order(tmp, decreasing=TRUE)]
    retval
  }

# general function
summary.genotype  <-  function(object,...,maxsum)
  {
    # if we are called from within summary.data.frame, fall back to
    # summary.factor so that we don't mess up the display
    if(!missing(maxsum))
      return(shortsummary.genotype(object,...,maxsum=maxsum))
    
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

