% groupGenotype.Rd
%--------------------------------------------------------------------------
% What: Group genotype values man page
% Time-stamp: <2007-07-21 12:02:09 ggorjan>
%--------------------------------------------------------------------------

\name{groupGenotype}

\alias{groupGenotype}

\title{Group genotype values}

\description{\code{groupGenotype} groups genotype or haplotype values
  according to given "grouping/mapping" information}

\usage{

groupGenotype(x, map, haplotype=FALSE, factor=TRUE, levels=NULL, verbose=FALSE)

}

\arguments{
  \item{x}{genotype or haplotype}
  \item{map}{list, mapping information, see details and examples}
  \item{haplotype}{logical, should values in a \code{map} be treated as
    haplotypes or genotypes, see details}
  \item{factor}{logical, should output be a factor or a character}
  \item{levels}{character, optional vector of level names if factor is
    produced (\code{factor=TRUE}); the default is to use the sort order
    of the group names in \code{map}}
  \item{verbose}{logical, print genotype names that match entries in
    the map - mainly used for debugging}
}

\value{A factor or character vector with genotypes grouped}

\details{

  Examples show how \code{map} can be constructed. This are the main
  points to be aware of:
  \itemize{
    \item names of list components are used as new group names

    \item list components hold genotype names per each group

    \item genotype names can be specified directly i.e. "A/B" or
    abbreviated such as "A/*" or even "*/*", where "*" matches any
    possible allele, but read also further on

    \item all genotype names that are not specified can be captured
    with ".else" (note the dot!)

    \item genotype names that were not specified (and ".else" was not
    used) are changed to \code{NA}
  }

  \code{map} is inspected before grouping of genotypes is being
  done. The following steps are done during inspection:
  \itemize{
    \item ".else" must be at the end (if not, it is moved) to match
    everything that has not yet been defined

    \item any specifications like "A/*", "*/A", or "*/*" are extended to
    all possible genotypes based on alleles in argument \code{alleles} -
    in case of \code{haplotype=FALSE}, "A/*" and "*/A" match the same
    genotypes

    \item since use of "*" and ".else" can cause duplicates along the
    whole map, duplicates are removed sequentially (first occurrence
    is kept)
  }

  Using ".else" or "*/*" at the end of the map produces the same result,
  due to removing duplicates sequentially.
}

\author{Gregor Gorjanc}

\seealso{
  \code{\link{genotype}},
  \code{\link{haplotype}},
  \code{\link{factor}}, and
  \code{\link{levels}}

}

\examples{

## --- Setup ---

x <- c("A/A", "A/B", "B/A", "A/C", "C/A", "A/D", "D/A",
       "B/B", "B/C", "C/B", "B/D", "D/B",
       "C/C", "C/D", "D/C",
       "D/D")
g <- genotype(x, reorder="yes")
## "A/A" "A/B" "A/B" "A/C" "A/C" "A/D" "A/D" "B/B" "B/C" "B/C" "B/D" "B/D"
## "C/C" "C/D" "C/D" "D/D"

h <- haplotype(x)
## "A/A" "A/B" "B/A" "A/C" "C/A" "A/D" "D/A" "B/B" "B/C" "C/B" "B/D" "D/B"
## "C/C" "C/D" "D/C" "D/D"

## --- Use of "A/A", "A/*" and ".else" ---

map <- list("homoG"=c("A/A", "B/B", "C/C", "D/D"),
            "heteroA*"=c("A/B", "A/C", "A/D"),
            "heteroB*"=c("B/*"),
            "heteroRest"=".else")

(tmpG <- groupGenotype(x=g, map=map, factor=FALSE))
(tmpH <- groupGenotype(x=h, map=map, factor=FALSE, haplotype=TRUE))

## Show difference between genotype and haplotype treatment
cbind(as.character(h), gen=tmpG, hap=tmpH, diff=!(tmpG == tmpH))
##              gen          hap          diff
##  [1,] "A/A" "homoG"      "homoG"      "FALSE"
##  [2,] "A/B" "heteroA*"   "heteroA*"   "FALSE"
##  [3,] "B/A" "heteroA*"   "heteroB*"   "TRUE"
##  [4,] "A/C" "heteroA*"   "heteroA*"   "FALSE"
##  [5,] "C/A" "heteroA*"   "heteroRest" "TRUE"
##  [6,] "A/D" "heteroA*"   "heteroA*"   "FALSE"
##  [7,] "D/A" "heteroA*"   "heteroRest" "TRUE"
##  [8,] "B/B" "homoG"      "homoG"      "FALSE"
##  [9,] "B/C" "heteroB*"   "heteroB*"   "FALSE"
## [10,] "C/B" "heteroB*"   "heteroRest" "TRUE"
## [11,] "B/D" "heteroB*"   "heteroB*"   "FALSE"
## [12,] "D/B" "heteroB*"   "heteroRest" "TRUE"
## [13,] "C/C" "homoG"      "homoG"      "FALSE"
## [14,] "C/D" "heteroRest" "heteroRest" "FALSE"
## [15,] "D/C" "heteroRest" "heteroRest" "FALSE"
## [16,] "D/D" "homoG"      "homoG"      "FALSE"

map <- list("withA"="A/*", "rest"=".else")
groupGenotype(x=g, map=map, factor=FALSE)
##  [1] "withA" "withA" "withA" "withA" "withA" "withA" "withA" "rest"  "rest"
## [10] "rest"  "rest"  "rest"  "rest"  "rest"  "rest"  "rest"

groupGenotype(x=h, map=map, factor=FALSE, haplotype=TRUE)
##  [1] "withA" "withA" "rest"  "withA" "rest"  "withA" "rest"  "rest"  "rest"
## [10] "rest"  "rest"  "rest"  "rest"  "rest"  "rest"  "rest"

## --- Use of "*/*" ---

map <- list("withA"="A/*", withB="*/*")
groupGenotype(x=g, map=map, factor=FALSE)
##  [1] "withA" "withA" "withA" "withA" "withA" "withA" "withA" "withB" "withB"
## [10] "withB" "withB" "withB" "withB" "withB" "withB" "withB"

## --- Missing genotype specifications produces NA's ---

map <- list("withA"="A/*", withB="B/*")
groupGenotype(x=g, map=map, factor=FALSE)
##  [1] "withA" "withA" "withA" "withA" "withA" "withA" "withA" "withB" "withB"
## [10] "withB" "withB" "withB" NA      NA      NA      NA

groupGenotype(x=h, map=map, factor=FALSE, haplotype=TRUE)
##  [1] "withA" "withA" "withB" "withA" NA      "withA" NA      "withB" "withB"
## [10] NA      "withB" NA      NA      NA      NA      NA

}

\keyword{manip}
\keyword{misc}

%--------------------------------------------------------------------------
% groupGenotype.Rd ends here