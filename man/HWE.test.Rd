% : HWE.test.Rd,v 1.13 2004/05/25 19:40:02 warnesgr Exp $

\name{HWE.test}
\alias{HWE.test}
\alias{HWE.test.genotype}
\alias{HWE.test.data.frame}
\alias{print.HWE.test}

\title{Estimate Disequilibrium and Test for Hardy-Weinberg Equilibrium}

\description{
  Estimate disequilibrium parameter and test the null hypothesis that 
  Hardy-Weinberg equilibrium holds. 
}

\usage{
HWE.test(x, ...)
\method{HWE.test}{genotype}(x, exact = nallele(x)==2, simulate.p.value=!exact,
         B=10000, conf=0.95, ci.B=1000, ... )
\method{HWE.test}{data.frame}(x, ..., do.Allele.Freq=TRUE, do.HWE.test=TRUE)
\method{print}{HWE.test}(x, show=c("D","D'","r","table"), ...)
}
\arguments{
  \item{x}{genotype or haplotype object.}
  \item{exact}{a logical value indicated whether the p-value should be
    computed using the exact method, which is only available for 2
    allele genotypes.}
  \item{simulate.p.value}{a logical value indicating whether the p-value
    should be computed using simulation instead of using the
    \eqn{\chi^2}{Chi-Square} approximation. Defaults to \code{TRUE}.}
  \item{B}{Number of simulation iterations to use when
    \code{simulate.p.value=TRUE}. Defaults to 10000.}
  \item{conf}{Confidence level to use when computing the confidence
    level for D-hat.  Defaults to 0.95, should be in (0,1). }
  \item{ci.B}{Number of bootstrap iterations to use when computing the
    confidence interval. Defaults to 1000.}
%   \item{ci.type}{Method of calculating the confidence interval using the
%     bootstrap sample.  Defaults to \code{"basic"}.  See 
%     \code{\link[boot]{boot.ci}} for details.}
  \item{show}{a character vector containing the names of HWE test
    statistics to display from the set of "D", "D'", "r", and "table".}
  \item{...}{ optional parameters passed to \code{HWE.test} (data.frame
    method) or \code{chisq.test} (base method).}
  \item{do.Allele.Freq}{logicial indication whether to summarize allele
    frequencies.}
  \item{do.HWE.test}{logicial indication whether to perform HWE tests}
}
\details{
  HWE.test calls \code{\link{diseq}} to computes the Hardy-Weinberg
  (dis)equilibrium statistics D, D', and r (correlation coefficient).
  Next it calls \code{\link{diseq.ci}} to compute a bootstrap confidence
  interval for these estimates.  Finally, it calls
  \code{\link{chisq.test}} to compute a p-value for Hardy-Weinberg
  Equilibrium using a simulation/permutation method.

  Using bootstrapping for the confidence interval and simulation for the
  p-value avoids reliance on the assumptions the underlying Chi-square
  approximation.  This is particularly important when some allele pairs
  have small counts.

  For details on the definition of D, D', and r, see the help page for
  \code{\link{diseq}}.
}
\value{
  An object of class \code{HWE.test} with components
  \item{diseq}{A \code{\link{diseq}} object providing details on the
    disequilibrium estimates.}
  \item{ci}{A \code{\link{diseq.ci}} object providing details on the
    bootstrap confidence intervals for the disequilibrium estimates.}
  \item{test}{A \code{htest} object providing details on the
    permutation based Chi-square test.}
  \item{call}{function call used to creat this object.}
  \item{conf, B, ci.B, simulate.p.value}{values used for these arguments.}
}
\author{ Gregory R. Warnes \email{greg@warnes.net }
  }

\seealso{
  \code{\link{genotype}},
  \code{\link{diseq}},
  \code{\link{diseq.ci}},
  \code{\link{HWE.chisq}},
  \code{\link{HWE.exact}},
  \code{\link[stats]{chisq.test}}
}


\examples{
\testonly{
set.seed(4657613)
}
example.data   <- c("D/D","D/I","D/D","I/I","D/D",
                    "D/D","D/D","D/D","I/I","")
g1  <- genotype(example.data)
g1

HWE.test(g1)

#compare with
diseq(g1)
diseq.ci(g1)
HWE.chisq(g1)
HWE.exact(g1)


three.data   <- c(rep("A/A",8),
                  rep("C/A",20),
                  rep("C/T",20),
                  rep("C/C",10),
                  rep("T/T",3))

g3  <- genotype(three.data)
g3

HWE.test(g3, ci.B=10000)


}
\keyword{ misc }
%%\keyword{genetics}
