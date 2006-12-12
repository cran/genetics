% $Id: power.casectl.Rd 173 2004-09-21 21:45:45Z warnes $

\name{power.casectrl}
\alias{power.casectrl}
\title{Power for case-control genetics study}
\description{
  Calculate power for case-control genetics study
}
\usage{
power.casectrl(N, gamma = 4.5, p = 0.15, kp = 0.1, alpha = 0.05, fc = 0.5, 
                minh = c("multiplicative", "dominant", "recessive"))
}
\arguments{
  \item{N}{total number of subjects}
  \item{gamma}{Relative risk in multiplicative model.  Not used in Dominant 
	or Recessive model.}
  \item{p}{frequency of A (protective) allele}
  \item{kp}{ significance level}
  \item{alpha}{ prevalence of disease }
  \item{fc}{ fraction of cases }
  \item{minh}{ mode of inheritance, one of  }
}
%\details{
%}
\value{
  power for the specified parameter values.
}
\references{Long, A. D. and C. H. Langley (1997). Genetic analysis of 
	    complex traits. Science 275: 1328.}
\author{Michael Man}
%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{

# single calc
power.casectrl(p=0.1, N=50, gamma=1.1, kp=.1, alpha=5e-2, minh='r')

# for a range of sample sizes
power.casectrl(p=0.1, N=c(25,50,100,200,500), gamma=1.1, kp=.1,
              alpha=5e-2, minh='r')

# create a power table
fun <- function(p)
  power.casectrl(p=p, N=seq(100,1000,by=100), gamma=1.1, kp=.1,
                alpha=5e-2, minh='recessive')

m <- sapply( X=seq(0.1,0.9, by=0.1), fun )
colnames(m) <- seq(0.1,0.9, by=0.1)
rownames(m) <- seq(100,1000,by=100)

print(round(m,2))
}
\keyword{design}
%\keyword{genetics}

