% $Id: HWE.exact.Rd,v 1.5 2004/09/21 21:45:45 warnes Exp $

\name{HWE.exact}
\alias{HWE.exact}
\title{Exact Test of Hardy-Weinberg Equilibrium for 2-Allele Markers}
\description{
  Exact test of Hardy-Weinberg Equilibrium for 2 Allele Markers.
}
\usage{
HWE.exact(x)
}
\arguments{
  \item{x}{ Genotype object }
}
%\details{
%}
\value{
  Object of class 'htest'.
}
\references{
    Emigh TH. (1980) "Comparison of tests for Hardy-Weinberg Equilibrium",
    Biometrics, 36, 627-642.
  }
\author{
  David Duffy \email{davidD@qimr.edu.au} with modifications by Gregory
  R. Warnes \email{gregory\_r\_warnes@groton.pfizer.com}
  }
\note{ This function only works for genotypes with exactly 2 alleles.}
\seealso{
  \code{\link{HWE.chisq}},
  \code{\link{HWE.test}},
  \code{\link{diseq}},
  \code{\link{diseq.ci}}
}

\examples{
example.data   <- c("D/D","D/I","D/D","I/I","D/D",
                    "D/D","D/D","D/D","I/I","")
g1  <- genotype(example.data)
g1

HWE.exact(g1)
# compare with
HWE.chisq(g1)


\testonly{
set.seed(465764)
}

g2 <- genotype(sample( c("A","C"), 100, p=c(100,10), rep=TRUE),
               sample( c("A","C"), 100, p=c(100,10), rep=TRUE) )
HWE.exact(g2)

}
\keyword{ misc }
%%\keyword{genetics}

