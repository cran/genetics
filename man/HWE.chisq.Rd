% $Id: HWE.chisq.Rd,v 1.4 2004/11/15 15:34:59 warnes Exp $

\name{HWE.chisq}
\alias{HWE.chisq}
\alias{HWE.chisq.genotype}

\title{Perform Chi-Square Test for Hardy-Weinberg Equilibrium}

\description{
  Test the null hypothesis that Hardy-Weinberg equilibrium holds using
  the Chi-Square method.
}

\usage{
HWE.chisq(x, ...)
\method{HWE.chisq}{genotype}(x, simulate.p.value=TRUE, B=10000, ...)
          
}
\arguments{
  \item{x}{genotype or haplotype object.}
  \item{simulate.p.value}{a logical value indicating whether the p-value
    should be computed using simulation instead of using the
    \eqn{\chi^2}{Chi-Square} approximation. Defaults to \code{TRUE}.}
  \item{B}{Number of simulation iterations to use when
    \code{simulate.p.value=TRUE}. Defaults to 10000.}
  \item{...}{ optional parameters passed to \code{chisq.test}}
}
\details{
  This function generates a 2-way table of allele counts, then calls
  \code{\link{chisq.test}} to compute a p-value for Hardy-Weinberg
  Equilibrium.  By default, it uses an unadjusted Chi-Square test
  statistic and computes the p-value using a simulation/permutation
  method.  When \code{simulate.p.value=FALSE}, it computes the test
  statistic using the Yates continuity correction and tests it against
  the asymptotic Chi-Square distribution with the approproate degrees of
  freedom.

  Note: The Yates continuty correction is applied *only* when
  \code{simulate.p.value=FALSE}, so that the reported test statistics
  when \code{simulate.p.value=FALSE} and \code{simulate.p.value=TRUE}
  will differ.
}
\value{
  An object of class \code{htest}.
}

\seealso{
  \code{\link{HWE.exact}},
  \code{\link{HWE.test}},
  \code{\link{diseq}},
  \code{\link{diseq.ci}},
  \code{\link{allele}},
  \code{\link{chisq.test}},
  \code{\link[boot]{boot}},
  \code{\link[boot]{bootci}}
}


\examples{
\testonly{
set.seed(4657613)
}
example.data   <- c("D/D","D/I","D/D","I/I","D/D",
                    "D/D","D/D","D/D","I/I","")
g1  <- genotype(example.data)
g1

HWE.chisq(g1)
# compare with
HWE.exact(g1)
# and 
HWE.test(g1)

three.data   <- c(rep("A/A",8),
                  rep("C/A",20),
                  rep("C/T",20),
                  rep("C/C",10),
                  rep("T/T",3))

g3  <- genotype(three.data)
g3

HWE.chisq(g3, B=10000)


}
\keyword{ misc }
%%\keyword{genetics}
