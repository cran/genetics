
#diseq.old <- function(x)
#{
#  if(!("genotype") %in% class(x) )
#    stop("x must inherit from class 'genotype'.")
  
#  # Estimates and tests per Wier (1996)  Genetic Data Analysis II

#  allele.names  <-  allele.names(x)
  
#  n <- length(na.omit(x))

#  require(ctest)

#  # specify levels so we get zeros where we need them
#  tab  <- table( factor(allele(x,1), levels=allele.names(x)),
#                 factor(allele(x,2), levels=allele.names(x)) )

#  allele.prob  <- table( allele(x) ) / (2 * n)

#  k  <-  length(allele.names)
#  D.hat  <-  matrix(NA, ncol=5, nrow= k * (k+1) / 2 + 1)
#  rnames  <- rep("",nrow(D.hat))
#  index  <- 1
#  for(i in 1:k)
#    for(j in i:k)
#      {
#        rnames[index]  <- paste( allele.names[i],
#                                allele.names[j], sep="/")

#        D.hat[index,1]  <- P  <-  tab[i,j]

#        if(i==j)
#          {
#            D.hat[index,2] <- pp <- n*allele.prob[allele.names[i]] ^ 2
#            D.hat[index,5] <- 1
#          }
#        else
#          {
#            D.hat[index,2] <- pp <- 2 * n * allele.prob[allele.names[i]] *
#                                            allele.prob[allele.names[j]]
#            D.hat[index,5] <- 2

#            D.hat[index,3] <- D  <- P - pp
#            D.hat[index,4] <- D/n / D.hat[index,5]
#          }
          
#        index  <-  index+1
#      }

#  rnames[nrow(D.hat)]  <- "Overall"
#  D.hat[nrow(D.hat),1]  <- n

#  ab  <- abs(D.hat[,3]/D.hat[,5])[-nrow(D.hat)]
##  ab  <- (D.hat[,3]/D.hat[,5])[-nrow(D.hat)]
#  nab  <- D.hat[,5][-nrow(D.hat)]
#  Dh  <- mean(rep(ab,nab),na.rm=TRUE)/n
#  D.hat[nrow(D.hat),4]  <-Dh

#  rownames(D.hat)  <- rnames
#  colnames(D.hat) <- c("Observed","Expected","Obs-Exp","D-hat","")
  
#  retval  <- list()
#  retval$data <-  tab
#  retval$D.hat  <- D.hat
#  retval$call  <- match.call()
#  class(retval)  <-  "diseq.old"
#  return(retval)
#}


#print.diseq.old  <-  function(x, show.table=TRUE, ...)
#  {

#    cat("\n")
#    if(!is.null(x$locus))
#      {
#        cat("\n")
#        print( x$locus )
#      }
#    cat("\n")
#    cat("Call: \n")
#    print(x$call)
#    cat("\n")
#    if(show.table)
#      {
#        cat("Disequlibrium (D-hat) Computation Table:\n")
#        cat("\n") 
#        print(x$D.hat[,-5] )
#        cat("\n")
#      }
#    cat("Overall Disequlibrium:\n")
#    cat("\n")
#    cat("\tD-hat :  ", x$D.hat[nrow(x$D.hat),4], "\n", sep="")
#    cat("\n")
#    cat("\n")
#  }

#diseq.ci.old <- function(x, R=1000, conf=0.95)
#{
#  if (!("genotype") %in% class(x) )
#    stop("x must inherit from class 'genotype'.")
  
#  bootfun <- function(x, ids) {
#    tmp <- diseq.old(x[ids])$D.hat
#    tmp[nrow(tmp),4]
#  }
  
#  bb <- boot( x, bootfun, R=R )
#  bb.ci <- boot.ci(bb, type=type, conf=conf, ...)
#  bb.ci
#}


diseq <- function(x, ...)
{
	UseMethod("diseq")
}

# Pairwise Disequilibrium Measures. For each pair of markers we
# calculated D, D' and r, the most commonly used measures of LD
# (other measures are reviewed in refs 14,15 ). For a pair of markers
# i and j, we defined Dij = p(11) - p(1·)p(·1), where p(ab) is the
# estimated frequency of the haplotypes with alleles a at marker i and
# b at marker j and · denotes any allele. Then, D'ij = |D11/Dmax|,
# where Dmax = max(p(1·)p(·1),p(2·)p(·2)) if Dij > 0 and Dmax =
# max(p(1·)p(·2),p(2·)p(·1)) otherwise, and r²ij = Dij²/(p(1·)p(·1)p(2·)p(·2)).

diseq.genotype <- function(x, ...)
  { 
    observed.no <- table( factor(allele(x,1), levels=allele.names(x)),
                          factor(allele(x,2), levels=allele.names(x)) )
    observed <- prop.table(observed.no)
    observed <- 1/2 * (observed + t(observed) )

    retval <- diseq.table(observed)
    retval$observed.no <- observed.no
    retval$call <- match.call()
    retval
  }
  
diseq.table <- function(x, ...)
{
  observed <- x
  allele.freq <- apply(observed,1,sum)
  # equal to: allele.freq <- apply(observed,2,sum)

  expected <- outer(allele.freq, allele.freq, "*")

  diseq <- expected - observed
  diag(diseq) <- NA

  dmax.positive <- expected
  # equals: max( p(i)p(j), p(j)p(i) )
  
  dmax.negative <- outer(allele.freq, allele.freq, pmin ) - expected
  # equals: min( p(i) * (1 - p(j)), p(j)( 1 - (1-p(i) ) ) )
  
  dprime <- diseq / ifelse( diseq > 0, dmax.positive, dmax.negative )

  # r gives the pairwise correlation coefficient for pairs containing at lease
  # one allele from the specified pair.
  # For two alleles:
  #    corr coefficient = diseq / sqrt( p(a) * (1-p(a) ) * p(b) * (1-p(b)) )
  p.1.minus.p <- allele.freq * (1-allele.freq)
  r <-  -diseq  / sqrt( outer( p.1.minus.p, p.1.minus.p, "*") )

  # above formula works unchanged for 2 alleles, but requires adjustment
  # for multiple alleles.  
  r <- r * (length(allele.freq) - 1)

  offdiag.expected <- expected
  diag(offdiag.expected) <- NA
  sum.expected <- sum(offdiag.expected, na.rm=TRUE)

  if(all(dim(x)==2)) # 2 allele case
    {
      diseq.overall <- diseq[1,2]
      dprime.overall <- dprime[1,2]
      r.overall <- r[1,2]
    }
  else
    {
      diseq.overall <- sum( abs(diseq) * expected , na.rm=TRUE ) / sum.expected
      dprime.overall <- sum( abs(dprime) * expected , na.rm=TRUE ) / sum.expected
      r.overall <- sum( abs(r) * expected , na.rm=TRUE ) / sum.expected
    }
  
  diag(r) <- 1.0

  retval <- list(
                 call = match.call(),
                 observed=observed,
                 expected=expected,
                 allele.freq=allele.freq,
                 D=diseq,
                 Dprime=dprime,
                 r=r,
                 D.overall=diseq.overall,
                 Dprime.overall=dprime.overall,
                 r.overall = r.overall
                 )

  class(retval) <- "diseq"
  retval
}


print.diseq  <-  function(x, show=c("D","D'","r","table"), ...)
  {

    cat("\n")
    if(!is.null(x$locus))
      {
        cat("\n")
        print( x$locus )
      }
    cat("\n")
    cat("Call: \n")
    print(x$call)
    cat("\n")
    if("D" %in% show)
      {
        cat("Disequlibrium for each allele pair (D)\n")
        cat("\n") 
        print(x$D)
        cat("\n")
      }
    if("D'" %in% show)
      {
        cat("Disequlibrium for each allele pair (D')\n")
        cat("\n") 
        print(x$Dprime)
        cat("\n")
      }
    if("r" %in% show)
      {
        cat("Correlation coefficient for each allele pair (r)\n")
        cat("\n") 
        print(x$r)
        cat("\n")
      }
    if("table" %in% show)
      {
        cat("Observed vs Expected frequency table\n")
        cat("\n")
        table <- cbind(
                       Obs=c(x$observed),
                       Exp=c(x$expected),
                       "Obs-Exp"= c(x$observed - x$expected)
                       )
        rownames(table) <- outer(rownames(x$observed),
                                 rownames(x$observed), paste, sep="/")
        
        print (table)
        cat("\n")
      }
    
    if( any(c("D","D'","r") %in% show))
      {
        if( ncol(x$r) <= 2 )
          cat("Overall Values\n")
        else
          cat("Overall Values (mean absolute-value weighted by expected allele frequency)\n")
        cat("\n")
        
        if("D" %in% show)
          cat("  D :  ", x$D.overall, "\n", sep="")
        if("D'" %in% show)
          cat("  D':  ", x$Dprime.overall, "\n", sep="")
        if("r" %in% show)
          cat("  r :  ", x$r.overall, "\n", sep="")
        cat("\n")
      }
    
    cat("\n")
  }

diseq.ci <- function(x, R=1000, conf=0.95, correct=TRUE, na.rm=TRUE, ...)
{
  if (!("genotype") %in% class(x) )
    stop("x must inherit from class 'genotype'.")

  if( any(is.na(x) ) )
    {
      if( na.rm)
        x <- na.omit(x)
      else
        stop("Missing values and NaN's not allowed if `na.rm' is FALSE.")
    }

  if( !require(combinat) )
    stop("Depends on availability of 'combinat' library")
  
  # step 1 - generate summary table
  observed.no <- table( factor(allele(x,1), levels=allele.names(x)),
                        factor(allele(x,2), levels=allele.names(x)) )
  observed <- prop.table(observed.no)
  observed <- 1/2 * (observed + t(observed) )

  # step 2 - make table into a probability vector for calling rmultinom
  n <- sum(observed.no)
  prob.vector <- c(observed)

  # step 3 - sample R multinomials with the specified frequenceis
  # (include observed data to avoid bias)
  resample.data <- cbind(c(observed.no),
                         rmultz2( n, prob.vector, R ) )
  
  bootfun <- function(x) {
    observed[,] <- x/n
    observed <- 1/2 * (observed + t(observed) )
    d <-  diseq(observed)
    c( "Overall D "=d$D.overall,
       "Overall D'"=d$Dprime.overall,
       "Overall r "=d$r.overall)
  }

  results <- apply( resample.data, 2, bootfun )

  alpha.2 <- (1-conf)/2

#  ci <- t(apply(results, 1,
#              quantile, c( alpha.2 , 1-alpha.2), na.rm=TRUE ))

  if(length(allele.names(x))<=2)
    {
      ci <- t(apply(results, 1, function(x) quantile(x, c(0.025, 0.975),
                                                     na.rm=na.rm ) ) )
      warning.text <- NULL
    }
  else
    {
      warning.text <- paste("For more than two alleles, overall",
                            "disequlibrium statistics are bounded",
                            "between [0,1].  Because of this, confidence",
                            "intervals for values near 0 and 1 are",
                            "ill-behaved.", sep=" ")
      
      if(correct)
        {
          warning.text <- paste(warning.text, "A rough correction has been applied, but",
                                "the intervals still may not be correct for values near 0 or 1.",
                                sep=" ")

          ci <- t(apply(results, 1,
                        function(x)
                        ci.balance(x,x[1],confidence=conf,
                                   minval=0,maxval=1)$ci ))
        }
      else
        ci <- t(apply(results, 1, function(x) quantile(x, c(0.025, 0.975) ) ) )
      
      warning(paste(strwrap(c(warning.text,"\n"),prefix="  "),collapse="\n") )
    }

  na.count <-  function(x) sum(is.na(x))
  nas <- apply( results, 1, na.count)

  zero.in.range <- (ci[,1] <= 0) & (ci[,2] >= 0)
  
  ci <- cbind( "Observed"=results[,1], ci, "NAs"=nas,
               "Zero in Range"=zero.in.range )

  outside.ci <- (ci[,1] < ci[,2]) | (ci[,1] > ci[,3])
  
  if( any(outside.ci) )
    warning("One or more observed value outide of confidence interval. Check results.")

  if(any(nas>0))
    warning("NAs returned from diseq call")


  retval <- list(
         call=match.call(),
         R=R,
         conf=conf,
         ci=ci,
         warning.text=warning.text
         )
}

