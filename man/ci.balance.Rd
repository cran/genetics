% $Id: ci.balance.Rd,v 1.7 2004/09/21 21:45:45 warnes Exp $

\name{ci.balance}
\alias{ci.balance}
\title{Experimental Function to Correct Confidence Intervals At or Near
  Boundaries of the Parameter Space by 'Sliding' the Interval on the
  Quantile Scale.}
\description{Experimental function to correct confidence intervals at or near
  boundaries of the parameter space by 'sliding' the interval on the
  quantile scale.}
\usage{
ci.balance(x, est, confidence=0.95, alpha=1-confidence, minval, maxval,
           na.rm=TRUE)
}
\arguments{
  \item{x}{Bootstrap parameter estimates.}
  \item{est}{Observed value of the parameter.}
  \item{confidence}{Confidence level for the interval. Defaults to
    0.95.}
  \item{alpha}{Type I error rate (size) for the interval. Defaults to
               1-\code{confidence}.}  
  \item{minval}{A numeric value specifying the lower bound of the
                parameter space. Leave unspecified (the default) if
                there is no lower bound.}
  \item{maxval}{A numeric value specifying the upper bound of the
                parameter space. Leave unspecified (the default) if
                there is no upper bound.}
  \item{na.rm}{ logical. Should missing values be removed?}
}
\details{
  EXPERIMENTAL FUNCTION:

  This function attempts to compute a proper \code{conf}*100\%
  confidence interval for parameters at or near the boundary of the
  parameter space using bootstrapped parameter estimates by 'sliding'
  the confidence interval on the quantile scale.
  
  This is accomplished by attempting to place a \code{conf} *100\%
  interval symmetrically *on the quantile scale* about the observed
  value.  If a symmetric interval would exceed the observed data at the
  upper (lower) end, a one-sided interval is computed with the upper
  (lower) boundary fixed at the the upper (lower) boundary of the
  parameter space.
}
\value{
  A list containing:
  \item{ci}{A 2-element vector containing the lower and upper confidence
    limits.  The names of the elements of the vector give the actual
    quantile values used for the interval or one of the character
    strings "Upper Boundary" or "Lower Boundary".}
  \item{overflow.upper, overflow.lower}{The number of elements beyond
    those observed that would be needed to compute a symmetric (on the
    quantile scale) confidence interval.}
  \item{n.above, n.below}{The number of bootstrap values which are above
    (below) the observed value.}
  \item{lower.n, upper.n}{The index of the value used for the endpoint
    of the confidence interval or the character string "Upper Boundary"
    ("Lower Boundary").}
}
\author{ Gregory R. Warnes \email{Gregory\_R\_Warnes@groton.pfizer.com }
}

\seealso{
  \code{\link[boot]{boot}},
  \code{\link[bootstrap]{bootstrap}},
  Used by \code{\link{diseq.ci}}.
}

\examples{
# These are nonsensical examples which simply exercise the
# computation. See the code to diseq.ci for a real example.
#
# FIXME: Add real example using boot or bootstrap.  

set.seed(7981357)
x <- abs(rnorm(100,1))
ci.balance(x,1, minval=0)
ci.balance(x,1)

x <- rnorm(100,1)
x <- ifelse(x>1, 1, x)
ci.balance(x,1, maxval=1)
ci.balance(x,1)
}
\keyword{misc}
%%\keyword{genetics}

