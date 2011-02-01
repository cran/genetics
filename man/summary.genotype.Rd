% $Id: summary.genotype.Rd 1346 2011-02-02 02:21:44Z warnes $

\name{summary.genotype}

\alias{summary.genotype}
\alias{print.summary.genotype}

\title{ Allele and Genotype Frequency from a Genotype or
  Haplotype Object}

\description{
    \code{summary.genotype} creates an object containing allele and
    genotype frequency from a \code{genotype} or \code{haplotype}
    object.  \code{print.summary.genotype} displays a
    \code{summary.genotype} object.
}

\usage{
  \method{summary}{genotype}(object, ..., maxsum)
  \method{print}{summary.genotype}(x,...,round=2)
}

\arguments{

  \item{object, x}{ an object of class \code{genotype} or \code{haplotype} (for
    \code{summary.genotype}) or an object of class
    \code{summary.genotype} (for \code{print.summary.genotype}) }
  \item{\dots}{ optional parameters.  Ignored by \code{summary.genotype},
    passed to \code{print.matrix} by \code{print.summary,genotype}.}
  \item{maxsum}{ specifying any value for the parameter
    maxsum will cause \code{summary.genotype} to fall back to
    \code{summary.factor}.}
  \item{round}{ number of digits to use when displaying proportions.}
 }

\details{
 Specifying any value for the parameter \code{maxsum} will cause fallback
 to \code{summary.factor}.  This is so that the function
 \code{summary.dataframe} will give reasonable output when it contains a
 genotype column.  (Hopefully we can figure out something better to do
 in this case.)
}


\value{
  The returned value of \code{summary.genotype} is an object of class
  \code{summary.genotype} which
  is a list with the following components:

  \item{locus }{locus information field (if present) from \code{x}}.

  \item{allele.names}{ vector of allele names }

  \item{allele.freq }{
    A two column matrix with one row for each allele, plus one row for
    \code{NA} values (if present).  The first column, \code{Count},
    contains the frequency of the corresponding allele value.  The
    second column, \code{Proportion}, contains the fraction of alleles
    with the corresponding allele value.  Note each observation contains
    two alleles, thus the \code{Count} field sums to twice the number of
    observations.
  }

  \item{genotype.freq}{
    A two column matrix with one row for each genotype, plus one row for
    \code{NA} values (if present). The first column, \code{Count}, contains the
    frequency of the corresponding genotype.  The second column,
    \code{Proportion}, contains the fraction of genotypes with the
    corresponding value.
  }

  \code{print.summary.genotype} silently returns the object \code{x}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Gregory R. Warnes \email{greg@warnes.net} }
%\note{ ~~further notes~~ }
\seealso{
  \code{\link{genotype}},
  \code{\link{HWE.test}},
  \code{\link{allele}},
  \code{\link{homozygote}},
  \code{\link{heterozygote}},
  \code{\link{carrier}},
%  \code{\link{summary.genotype}},
  \code{\link{allele.count}}
  \code{\link{locus}}
  \code{\link{gene}}
  \code{\link{marker}}
  }

\examples{

example.data   <- c("D/D","D/I","D/D","I/I","D/D",
                    "D/D","D/D","D/D","I/I","")
g1  <- genotype(example.data)
g1

summary(g1)
}
\keyword{ misc }
%%\keyword{genetics}
