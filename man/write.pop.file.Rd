% $Id: write.pop.file.Rd,v 1.5 2004/09/17 18:46:26 warnes Exp $

\name{write.pop.file}
\alias{write.pop.file}
\alias{write.pedigree.file}
\alias{write.marker.file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create genetics data files}
\description{
  \code{write.pop.file} creates a 'pop' data file, as used by the
  GenePop (\url{http://wbiomed.curtin.edu.au/genepop/}) and LinkDos
  (\url{http://wbiomed.curtin.edu.au/genepop/linkdos.html}) software
  packages.

  \code{write.pedigree.file} creates a 'pedigree' data file, as used
  by the QTDT software package
  (\url{http://www.sph.umich.edu/statgen/abecasis/QTDT/}).
  
  \code{write.marker.file} creates a 'marker' data file, as used by
  the QTDT software package 
  (\url{http://www.sph.umich.edu/statgen/abecasis/QTDT/}).
}
\usage{
write.pop.file(data, file = "", digits = 2, description = "Data from R")
write.pedigree.file(data, family, pid, father, mother, sex,
                    file="pedigree.txt")
write.marker.file(data, location, file="marker.txt")
}
\arguments{
  \item{data}{Data frame containing genotype objects to be exported}
  \item{file}{Output filename}
  \item{digits}{Number of digits to use in numbering genotypes, either 2
    or 3.}
  \item{description}{Description to use as the first line of the 'pop'
    file.}
  \item{family, pid, father, mother}{Vector of family, individual,
    father, and mother id's, respectively.}
  \item{sex}{Vector giving the sex of the individual (1=Make, 2=Female)}
  \item{location}{Location of the marker relative to the gene of
    interest, in base pairs.}
}
\details{
  The format of 'Pop' files is documented at
  \url{http://wbiomed.curtin.edu.au/genepop/help_input.html}, the format
  of 'pedigree' files is documented at \url{http://www.sph.umich.edu/csg/abecasis/GOLD/docs/pedigree.html} and the format of 'marker'
  files is documented at
  \url{http://www.sph.umich.edu/csg/abecasis/GOLD/docs/map.html}.
}
\value{
  No return value.
}
\author{Gregory R. Warnes \email{gregory\_r\_warnes@groton.pfizer.com}}
\seealso{\code{\link{write.table}}}
\examples{
  # TBA
}
\keyword{IO}
