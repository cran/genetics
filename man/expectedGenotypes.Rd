\name{expectedGenotypes}
\alias{expectedGenotypes}
\title{Construct expected genotypes according to known allele variants}

\description{
  \code{expectedGenotypes} constructs expected genotypes according to 
  known allele variants, which can be quite tedious with large number of 
  allele variants. It can handle data with different level of ploidy.
}

\usage{
expectedGenotypes(x, alleles = allele.names(x), ploidy = 2)
}

\arguments{
  \item{x}{genotype object, as genotype.}
  \item{alleles}{vector of allele names, as character.}
  \item{ploidy}{number of chromosome sets i.e. 2 for human autosomal
                genes, as integer.}
  At least one of \code{x} or \code{alleles} must be given.
}

\value{
  A character vector with genotype names as "alele1/alele2" for diploid
  example.
}

\author{ Gregor GORJANC }

\seealso{ \code{\link{allele.names}}, \code{\link{genotype}}}

\examples{
  \dontrun{ Scrapie example }

  scrapie <- c("ARQ/ARQ", "ARQ/ARQ", "ARR/ARQ", "AHQ/ARQ", "ARQ/ARQ")
  expectedGenotypes(as.genotype(scrapie))
  expectedGenotypes(alleles=c("ARR", "AHQ", "ARH", "ARQ", "VRR", "VRQ"))  
  scrapie <- genotype(scrapie,
                      alleles=c("ARR", "AHQ", "ARH", "ARQ", "VRR", "VRQ"),
                      reorder="yes")
  expectedGenotypes(scrapie)
  
}

\keyword{ manip }
