% $Id: power.casectl.Rd 1272 2007-07-18 11:57:28Z ggorjan $

\name{power.casectrl}
\alias{power.casectrl}
\title{Power for case-control genetics study}
\description{
  Calculate power for case-control genetics study
}
\usage{
power.casectrl(N, gamma = 4.5, p = 0.15, kp = 0.1, alpha = 0.05, fc = 0.5,
                minh = c("multiplicative", "dominant", "recessive", "partialrecessive") )
}
\arguments{
  \item{N}{total number of subjects}
  \item{gamma}{Relative risk in multiplicative model.  Not used in Dominant
	or Recessive model.}
  \item{p}{frequency of A (protective) allele}
  \item{kp}{ significance level}
  \item{alpha}{ prevalence of disease }
  \item{fc}{ fraction of cases }
  \item{minh}{ mode of inheritance, one of "multiplicative", "dominant",
     "recessive", or "partialrecessive" }
}
\details{
\preformatted{
under HWE                       AA       Aa       aa
  fHW = p(genotype)        = c( p^2,     2pq,     q^2 )

model specification
  f.mod = relative risk    = c(gamma^2,  gamma,    1 )    #   multiplicative model
  f.mod =                  = c(  0,        0,      1 )    #         dominant model
  f.mod =                  = c(  0,        1,      1 )    #        recessive model
  f.mod =                  = c(gamma,      1,      1 )    # partialrecessive model
conditional prob.
  p(D|genotype) = f.mod*pi = c(gamma^2,  gamma,    1 )*pi

population joint prob. (f.mod = 1 under Ho)
  Case     p(D,     genotype) = p(genotype)*     p(D|genotype)  = fHW*   f.mod*pi
  Control  p(D_not, genotype) = p(genotype)*(1 - p(D|genotype)) = fHW*(1-f.mod*pi)

population conditional prob. (f.mod = 1 under Ho)
  Case     p(genotype|D)     = p(D    , genotype)/P(D    ) = P(D    , genotype)/sum(P(D    , genotype)) = fHW*   f.mod*pi  /    sum(fHW*f.mod*pi)
  Control  p(genotype|D_not) = p(D_not, genotype)/P(D_not) = P(D_not, genotype)/sum(P(D_not, genotype)) = fHW*(1-f.mod*pi) / (1-sum(fHW*f.mod*pi))

sample or allocation probability
  1:1 case-control design  p(D|Sample) = fc = 1/2
  1:2 case-control design                fc = 1/3
  a prospective design                   fc = sum(fHW*f.mod*pi)

sample joint prob. (f.mod = 1 under Ho)
for prospective design, this is the same as population joint prob. since 'fc' cancels out with 'sum(fHW*f.mod*pi)'
  Case     p(genotype,D    |sample) = p(genotype|D    )*     p(D|Sample)  =    fc *fHW*   f.mod*pi  /    sum(fHW*f.mod*pi)
  Control  p(genotype,D_not|sample) = p(genotype|D_not)*(1 - p(D|Sample)) = (1-fc)*fHW*(1-f.mod*pi) / (1-sum(fHW*f.mod*pi))
}

}
\value{
  power for the specified parameter values.
}
\references{Long, A. D. and C. H. Langley (1997). Genetic analysis of
	    complex traits. Science 275: 1328.}
\author{Michael Man /email{Michael.Man@pfizer.com}}
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