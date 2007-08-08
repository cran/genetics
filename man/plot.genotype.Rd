\name{plot.genotype}

\alias{plot.genotype}

\title{Plot genotype object}

\description{
  \code{plot.genotype} can plot genotype or allele frequency of a genotype
  object.
}

\usage{
\method{plot}{genotype}(x, type=c("genotype", "allele"),
  what=c("percentage", "number"), ...)
}

\arguments{
  \item{x}{genotype object, as genotype.}
  \item{type}{plot "genotype" or "allele" frequency, as character.}
  \item{what}{show "percentage" or "number", as character}
  \item{\dots}{Optional arguments for \code{barplot}.}
}

\value{
  The same as in \code{barplot}.
}

\author{ Gregor Gorjanc }

\seealso{ \code{\link{genotype}}, \code{\link{barplot}} }

\examples{
  set <- c("A/A", "A/B", "A/B", "B/B", "B/B", "B/B",
           "B/B", "B/C", "C/C", "C/C")
  set <- genotype(set, alleles=c("A", "B", "C"), reorder="yes")
  plot(set)
  plot(set, type="allele", what="number")
}

\keyword{hplot}
