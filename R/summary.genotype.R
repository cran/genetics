# $Id: summary.genotype.R,v 1.10 2003/08/04 13:48:40 warnesgr Exp $

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
    paf <- prop.table(af)
    retval$allele.freq    <- cbind("Count"=af,"Proportion"=paf)

    gf  <- table( object )
    pgf <- prop.table(gf)
    retval$genotype.freq  <- cbind("Count"=gf,"Proportion"=pgf)


    ### from code submitted by David Duffy <davidD@qimr.edu.au>
    #
    n.typed<-sum(gf)
    correction<-n.typed/max(1,n.typed-1)
    ehet<-(1-sum(paf*paf))
    matings<- (paf %*% t(paf))^2
    uninf.mating.freq <- sum(matings)-sum(diag(matings))
    pic<- ehet - uninf.mating.freq

    retval$Hu <- correction * ehet
    retval$pic <- pic
    retval$n.typed <- n.typed
    retval$n.total <- length(object)
    retval$nallele <- nallele(object)
    #
    ###

    
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
    cat("Number persons typed: ", x$n.typed,
      " (", round(100*x$n.typed/x$n.total,1), "%)\n", sep="")

    cat("\n")
    cat("Allele Frequency: (", x$nallele, " alleles)\n", sep="")
    print(round(x$allele.freq,digits=round),...)
    cat("\n")


    cat("\n")
    cat("Genotype Frequency:\n")
    print(round(x$genotype.freq,digits=round),...)
    cat("\n")
    
    cat("Heterozygosity (Hu)  = ",  x$Hu, "\n", sep="")
    cat("Poly. Inf. Content   = ",  x$pic, "\n", sep="")
    cat("\n")
    
    invisible(x)
  }

