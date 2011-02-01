% $Id: gregorius.Rd 1346 2011-02-02 02:21:44Z warnes $

\name{gregorius}
\alias{gregorius}
\title{Probability of Observing All Alleles with a Given Frequency in a
  Sample of a Specified Size.}
\description{
  Probability of observing all alleles with a given frequency in a
  sample of a specified size.
}
\usage{
gregorius(freq, N, missprob, tol = 1e-10, maxN = 10000, maxiter=100, showiter = FALSE)
}
\arguments{
  \item{freq}{(Minimum) Allele frequency (required)}
  \item{N}{Number of sampled genotypes}
  \item{missprob}{Desired maximum probability of failing to observe an allele.}
  \item{tol}{Omit computation for terms which contribute less than this value.}
  \item{maxN}{Largest value to consider when searching for N.}
  \item{maxiter}{Maximum number of iterations to use when searching for N.}
  \item{showiter}{Boolean flag indicating whether to show the iterations
    performed when searching for N.}
}
\details{
  If \code{freq} and \code{N} are provided, but \code{missprob} is omitted,
  this function computes the probability of failing to observe all alleles
  with true underlying frequency \code{freq} when \code{N} diploid
  genotypes are sampled.  This is accomplished using the sum provided in
  Corollary 2 of Gregorius (1980), omitting terms which contribute less
  than \code{tol} to the result.

  When \code{freq} and \code{missprob} are provide, but \code{N} is
  omitted. A binary search on the range of [1,\code{maxN}] is performed
  to locate the smallest sample size, \code{N}, for which the
  probability of failing to observe all alleles with true
  underlying frequency \code{freq} is at most \code{missprob}.  In this
  case, \code{maxiter} specifies the largest number of iterations to use
  in the binary search, and \code{showiter} controls whether the
  iterations of the search are displayed.
}
\value{
  A list containing the following values:
  \item{call}{   Function call used to generate this object.}
  \item{method}{ One of the strings, "Compute missprob given N and freq",
    or "Determine minimal N given missprob and freq", indicating which
    type of computation was performed.}
  \item{retval$freq}{ Specified allele frequency.}
  \item{retval$N}{    Specified or computed sample size. }
  \item{retval$missprob}{ Computed probability of failing to observe all
  of the alleles with frequency \code{freq}. }
}
\references{
  Gregorius, H.R. 1980. The probability of losing an allele when
  diploid genotypes are sampled.  Biometrics 36, 643-652.
}
\note{
  This code produces sample sizes that are slightly larger than those
  given in table 1 of Gregorius (1980).  This appears to be due to
  rounding of the computed \code{missprob}s by the authors of that
  paper.
  }
    
\author{ Code submitted by David Duffy \email{davidD@qumr.edu.au},
  substantially enhanced by Gregory R. Warnes
  \email{greg@warnes.net}. }
%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }
\examples{

# Compute the probability of missing an allele with frequency 0.15 when
# 20 genotypes are sampled:
gregorius(freq=0.15, N=20)

# Determine what sample size is required to observe all alleles with true
# frequency 0.15 with probability 0.95
gregorius(freq=0.15, missprob=1-0.95)


}
\keyword{misc}
