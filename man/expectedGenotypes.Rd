\name{expectedGenotypes}

\alias{expectedGenotypes}
\alias{expectedHaplotypes}

\title{Construct expected genotypes/haplotypes according to known allele
  variants}

\description{
  \code{expectedGenotypes} constructs expected genotypes according to
  known allele variants, which can be quite tedious with large number of
  allele variants. It can handle different level of ploidy.
}

\usage{
expectedGenotypes(x, alleles=allele.names(x), ploidy=2, sort=TRUE,
                  haplotype=FALSE)
expectedHaplotypes(x, alleles=allele.names(x), ploidy=2, sort=TRUE,
                   haplotype=TRUE)
}

\arguments{
  \item{x}{genotype or haplotype}
  \item{alleles}{character, vector of allele names}
  \item{ploidy}{numeric, number of chromosome sets i.e. 2 for human autosomal
    genes}
  \item{sort}{logical, sort genotypes according to order of alleles in
    \code{alleles} argument}
  \item{haplotype}{logical, construct haplotypes i.e. ordered genotype}
  At least one of \code{x} or \code{alleles} must be given.
}

\details{

\code{expectedHaplotypes()} just calls \code{expectedGenotypes()} with
argument \code{haplotype=TRUE}.

}

\value{

A character vector with genotype names as "alele1/alele2" for diploid
example. Length of output is \eqn{(n*(n+1))/2} for genotype (unordered
genotype) and \eqn{n*n} for haplotype (ordered genotype) for \eqn{n}
allele variants.

}

\author{ Gregor Gorjanc }

\seealso{ \code{\link{allele.names}}, \code{\link{genotype}}}

\examples{
  ## On genotype
  prp <- c("ARQ/ARQ", "ARQ/ARQ", "ARR/ARQ", "AHQ/ARQ", "ARQ/ARQ")
  alleles <- c("ARR", "AHQ", "ARH", "ARQ", "VRR", "VRQ")
  expectedGenotypes(as.genotype(prp))
  expectedGenotypes(as.genotype(prp, alleles=alleles))
  expectedGenotypes(as.genotype(prp, alleles=alleles, reorder="yes"))

  ## Only allele names
  expectedGenotypes(alleles=alleles)
  expectedGenotypes(alleles=alleles, ploidy=4)

  ## Haplotype
  expectedHaplotypes(alleles=alleles)
  expectedHaplotypes(alleles=alleles, ploidy=4)[1:20]
}

\keyword{ manip }
