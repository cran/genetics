# $Id: ci.balance.R,v 1.1 2002/11/08 19:53:57 warnesgr Exp $
#
# $Log: ci.balance.R,v $
# Revision 1.1  2002/11/08 19:53:57  warnesgr
#
# - Moved ci.balance() to a separate file and created a documenation file for it.
# - Modified ci.balance to better annotate when it uses boundary values.
# - Modified diseq.ci to better provide warning message when the number of
#   alleles is greater than 3.
#
#

ci.balance <- function(x, est, conf=0.95, minval, maxval, na.rm=TRUE)
  {
    if( any(is.na(x) ) )
      {
        if( na.rm)
          x <- na.omit(x)
        else
          stop("Missing values and NaN's not allowed if `na.rm' is FALSE.")
      }

    if(missing(minval))
      {
        minval <- min(x)
        minname <- "min(x)"
      }
    else
      minname <- "Lower Boundary"
    
    if(missing(maxval))
      {
        maxval <- max(x)
        maxname <- "max(x)"
      }
    else
      maxname <- "Upper Boundary"
    
    x <- sort(x)
    n <- length(x)
    half.window <- n * conf / 2
    n.below <- sum( x < est ) + sum( x==est )/2
    n.above <- sum( x > est ) + sum( x==est )/2 

    overflow.upper <- max(0, half.window - n.above )
    overflow.lower <- max(0, half.window - n.below ) 
    
    lower.n <- max(1, floor  ( n.below - half.window - overflow.upper ) )
    upper.n <- min(n, ceiling( n - (n.above - half.window - overflow.lower ) ) )
    
    ci <- c( x[lower.n], x[upper.n] )
    names(ci) <- paste( format( c(lower.n, upper.n)/n*100,digits=3 ), "%", sep="")

    if(overflow.lower>0)
      {
        lower.n <- minname
        names(ci)[1] <- minname
        ci[1] <- minval
      }
    if(overflow.upper>0)
      {
        upper.n <- maxname
        names(ci)[2] <- maxname
        ci[2] <- maxval
      }
    

    return(
           ci=ci,
           overflow.upper=overflow.upper,
           overflow.lower=overflow.lower,
           n.above=n.above,
           n.below=n.below,
           lower.n=lower.n,
           upper.n=upper.n
           )
  }
