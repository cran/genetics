% order.genotype.Rd
%--------------------------------------------------------------------------
% What: Order/sort genotype man page
% Time-stamp: <2007-07-20 03:11:39 ggorjan>
%--------------------------------------------------------------------------

\name{order.genotype}

\alias{order}
\alias{order.genotype}
\alias{sort.genotype}
\alias{genotypeOrder}
\alias{genotypeOrder<-}

\concept{order.haplotype}
\concept{sort.haplotype}

\title{Order/sort genotype/haplotype object}

\description{Order/sort genotype or haplotype object according
  to order of allele names or genotypes}

\usage{

\method{order}{genotype}(..., na.last=TRUE, decreasing=FALSE,
  alleleOrder=allele.names(x), genotypeOrder=NULL)

\method{sort}{genotype}(x, decreasing=FALSE, na.last=NA, ...,
  alleleOrder=allele.names(x), genotypeOrder=NULL)

genotypeOrder(x)
genotypeOrder(x) <- value

}

\arguments{
  \item{\ldots}{genotype or haplotype in \code{order} method; not used
    for \code{sort} method}
  \item{x}{genotype or haplotype in \code{sort} method}
  \item{na.last}{as in default \code{\link{order}} or \code{\link{sort}}}
  \item{decreasing}{as in default \code{\link{order}} or \code{\link{sort}}}
  \item{alleleOrder}{character, vector of allele names in wanted order}
  \item{genotypeOrder}{character, vector of genotype/haplotype names in
    wanted order}
  \item{value}{the same as in argument \code{order.genotype}}
}

\value{The same as in \code{order} or \code{sort}}

\details{

Argument \code{genotypeOrder} can be usefull, when you want that some
genotypes appear "together", whereas they are not "together" by allele
order.

Both methods (\code{order} and \code{sort}) work with genotype and
haplotype classes.

If \code{alleleOrder} is given, \code{genotypeOrder} has no effect.

Genotypes/haplotypes, with missing alleles in \code{alleleOrder} are
treated as \code{NA} and ordered according to \code{\link{order}}
arguments related to \code{NA} values. In such cases a warning is issued
("Found data values not matching specified alleles. Converting to NA.")
and can be safely ignored. Genotypes present in \code{x}, but not
specified in \code{genotypeOrder}, are also treated as \code{NA}.

Value of \code{genotypeOrder} such as "B/A" matches also "A/B" in case
of genotypes.

Only unique values in argument \code{alleleOrder} or
\code{genotypeOrder} are used i.e. first occurrence prevails.

}

\author{Gregor Gorjanc}

\seealso{
  \code{\link{genotype}},
  \code{\link{allele.names}},
  \code{\link{order}}, and
  \code{\link{sort}}

}

\examples{
  x <- c("C/C", "A/C", "A/A", NA, "C/B", "B/A", "B/B", "B/C", "A/C")
  alleles <- c("A", "B", "C")

  g <- genotype(x, alleles=alleles, reorder="yes")
  ## "C/C" "A/C" "A/A" NA    "B/C" "A/B" "B/B" "B/C" "A/C"

  h <- haplotype(x, alleles=alleles)
  ## "C/C" "A/C" "A/A" NA    "C/B" "B/A" "B/B" "B/C" "A/C"

  ## --- Standard usage ---

  sort(g)
  ## "A/A" "A/B" "A/C" "A/C" "B/B" "B/C" "B/C" "C/C" NA

  sort(h)
  ## "A/A" "A/C" "A/C" "B/A" "B/B" "B/C" "C/B" "C/C" NA

  ## --- Reversed order of alleles ---

  sort(g, alleleOrder=c("B", "C", "A"))
  ## "B/B" "B/C" "B/C" "A/B" "C/C" "A/C" "A/C" "A/A" NA
  ## note that A/B comes after B/C since it is treated as B/A;
  ## order of alleles (not in alleleOrder!) does not matter for a genotype

  sort(h, alleleOrder=c("B", "C", "A"))
  ## "B/B" "B/C" "B/A" "C/B" "C/C" "A/C" "A/C" "A/A" NA

  ## --- Missing allele(s) in alleleOrder ---

  sort(g, alleleOrder=c("B", "C"))
  ## "B/B" "B/C" "B/C" "C/C" "A/C" "A/A" NA    "A/B" "A/C"

  sort(g, alleleOrder=c("B"))
  ## "B/B" "C/C" "A/C" "A/A" NA    "B/C" "A/B" "B/C" "A/C"
  ## genotypes with missing allele are treated as NA

  sort(h, alleleOrder=c("B", "C"))
  ## "B/B" "B/C" "C/B" "C/C" "A/C" "A/A" NA    "B/A" "A/C"

  sort(h, alleleOrder=c("B"))
  ## "B/B" "C/C" "A/C" "A/A" NA    "C/B" "B/A" "B/C" "A/C"

  ## --- Use of genotypeOrder ---

  sort(g, genotypeOrder=c("A/A", "C/C", "B/B", "A/B", "A/C", "B/C"))
  ## "A/A" "C/C" "B/B" "A/B" "A/C" "A/C" "B/C" "B/C" NA

  sort(h, genotypeOrder=c("A/A", "C/C", "B/B",
                          "A/C", "C/B", "B/A", "B/C"))
  ## "A/A" "C/C" "B/B" "A/C" "A/C" "C/B" "B/A" "B/C" NA

  ## --- Missing genotype(s) in genotypeOrder ---

  sort(g, genotypeOrder=c(       "C/C",        "A/B", "A/C", "B/C"))
  ## "C/C" "A/B" "A/C" "A/C" "B/C" "B/C" "A/A" NA    "B/B"

  sort(h, genotypeOrder=c(       "C/C",        "A/B", "A/C", "B/C"))
  ## "C/C" "A/C" "A/C" "B/C" "A/A" NA    "C/B" "B/A" "B/B"
}

\keyword{manip}
\keyword{misc}

%--------------------------------------------------------------------------
% order.genotype.Rd ends here