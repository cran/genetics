# $Id: binsearch.R,v 1.2 2002/12/11 21:29:55 warnesgr Exp $
#
# $Log: binsearch.R,v $
# Revision 1.2  2002/12/11 21:29:55  warnesgr
# - fixed typo
#
# Revision 1.1  2002/12/11 21:07:37  warnesgr
#
# - Moved binsearch() into a separate file and added documentation.
# - Simplified and documented binsearch() code.
# - binsearch() now detects whether the function is increasing or
#   decreasing and acts accordingly.
# - Updated gregorius() to use the modified binsearch.
#
#


binsearch <- function(fun, range, ..., target=0,
                      lower=ceiling(min(range)),upper=floor(max(range)),
                      maxiter=100, showiter=FALSE)
    {

      # initialize
      lo <- lower
      hi <- upper
      counter <- 0
      val.lo <- fun(lo,...)
      val.hi <- fun(hi,...)

      # check whether function is increasing or decreasing, & set sign
      # appropriately.
      if( val.lo > val.hi )
        sign <- -1
      else
        sign <- 1

      # check if value is outside specified range
      if(target * sign < val.lo * sign)
          outside.range <- TRUE
      else if(target * sign >  val.hi * sign)
          outside.range <- TRUE
      else
        outside.range <- FALSE

      # iteratively move lo & high closer together until we run out of
      # iterations, or they are adjacent, or they are identical
      while(counter < maxiter && !outside.range )
        {

          counter <- counter+1

          if(hi-lo<=1 || lo<lower || hi>upper) break;

          center <- round((hi - lo)/2 + lo ,0  )
          val <- fun(center)

          if(showiter)
            {
              cat("--------------\n")
              cat("Iteration #", counter, "\n")
              cat("lo=",lo,"\n")
              cat("hi=",hi,"\n")
              cat("center=",center,"\n")
              cat("fun(lo)=",val.lo,"\n")
              cat("fun(hi)=",val.hi,"\n")
              cat("fun(center)=",val,"\n")
            }

          
          if( val==target )
            {
              val.lo <- val.hi <- val
              lo <- hi <- center
              break;
            }
          else if( sign*val < sign*target )
            {
              lo <- center
              val.lo <- val
            }
          else #( val > target )
            {
              hi <- center
              val.hi <- val
            }

        if(showiter)
          {
            cat("new lo=",lo,"\n")
            cat("new hi=",hi,"\n")
            cat("--------------\n")
          }
          
        }
      
      # Create return value
      retval <- list()
      retval$call <- match.call()
      retval$numiter <- counter

      if( outside.range )
        {
          if(target * sign < val.lo * sign)
            {
              warning("Reached lower boundary")
              retval$flag="Lower Boundary"
              retval$where=lo
              retval$value=val.lo
            }
          else #(target * sign >  val.hi * sign)
          {
            warning("Reached upper boundary")
            retval$flag="Upper Boundary"
            retval$where=hi
            retval$value=val.hi
          }
        }
      else if( counter >= maxiter )
        {
          warning("Maximum number of iterations reached")
          retval$flag="Maximum number of iterations reached"
          retval$where=c(lo,hi)
          retval$value=c(fun.lo,fun.hi)
        }
      else if( val.lo==target )
        {
          retval$flag="Found"
          retval$where=lo
          retval$value=val.lo
        }
      else if( val.hi==target )
        {
          retval$flag="Found"
          retval$where=lo
          retval$value=val.lo
        }
      else
        {
          retval$flag="Between Elements"
          retval$where=c(lo, hi)
          retval$value=c(val.lo, val.hi)
        }
      
      return(retval)

    }
                 
                                             

