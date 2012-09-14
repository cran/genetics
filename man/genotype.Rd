% $Id: genotype.Rd 1353 2012-08-14 14:23:17Z warnes $

\name{genotype}
\alias{genotype}
\alias{haplotype}
\alias{is.genotype}
\alias{is.haplotype}
\alias{as.genotype}
\alias{as.haplotype}

%%
\alias{print.genotype}
\alias{==.genotype}
\alias{==.haplotype}
\alias{\%in\%}
\alias{\%in\%.default}
\alias{\%in\%.genotype}
\alias{\%in\%.haplotype}
\alias{[.genotype}
\alias{[.haplotype}
\alias{[<-.genotype}
\alias{[<-.haplotype}
\alias{heterozygote.genotype}
\alias{homozygote.genotype}
\alias{print.allele.count}
\alias{print.allele.genotype}
\alias{allele.count.genotype}
%\alias{allele.genotype}
%\alias{allele.names}
\alias{as.genotype.allele.count}
\alias{as.genotype.character}
\alias{as.genotype.default}
\alias{as.genotype.factor}
\alias{as.genotype.genotype}
\alias{as.genotype.haplotype}
\alias{as.genotype.table}
\alias{nallele}

\title{Genotype or Haplotype Objects.}
\description{
  \code{genotype} creates a genotype object.

  \code{haplotype} creates a haplotype object.

  \code{is.genotype} returns \code{TRUE} if \code{x} is of class
  \code{genotype}

  \code{is.haplotype} returns \code{TRUE} if \code{x} is of class
  \code{haplotype}

  \code{as.genotype} attempts to coerce its argument into an object of
  class \code{genotype}.

  \code{as.genotype.allele.count} converts allele counts (0,1,2) into
  genotype pairs ("A/A", "A/B", "B/B").

  \code{as.haplotype} attempts to coerce its argument into an object of
  class \code{haplotype}.

  \code{nallele} returns the number of alleles in an object of class
  \code{genotype}.

}
\usage{
  genotype(a1, a2=NULL, alleles=NULL, sep="/", remove.spaces=TRUE,
           reorder = c("yes", "no", "default", "ascii", "freq"),
           allow.partial.missing=FALSE, locus=NULL,
           genotypeOrder=NULL)

  haplotype(a1, a2=NULL, alleles=NULL, sep="/", remove.spaces=TRUE,
            reorder="no", allow.partial.missing=FALSE, locus=NULL,
            genotypeOrder=NULL)

  is.genotype(x)

  is.haplotype(x)

  as.genotype(x, ...)

  \method{as.genotype}{allele.count}(x, alleles=c("A","B"), \dots )

  as.haplotype(x, ...)

  \method{print}{genotype}(x, \dots)

  nallele(x)
}

\arguments{
  \item{x}{ either an object of class \code{genotype} or
    \code{haplotype} or an object to be converted to class \code{genotype} or
    \code{haplotype}.}
  \item{a1,a2}{ vector(s) or matrix containing two alleles
    for each individual. See details, below.}
  \item{alleles}{ names (and order if \code{reorder="yes"}) of possible
    alleles.}
  \item{sep}{ character separator or column number used to divide
    alleles when \code{a1} is a vector of strings where each string
    holds both alleles. See below for details.}
  \item{remove.spaces}{ logical indicating whether spaces and tabs will
    be removed from a1 and a2  before processing.}
  \item{reorder}{how should alleles within an individual be reordered.
    If \code{reorder="no"}, use the order specified by the alleles
    parameter.  If \code{reorder="freq"} or
    \code{reorder="yes"}, sort alleles within each individual by
    observed frequency.  If \code{reorder="ascii"}, reorder alleles in
    ASCII order (alphabetical, with all upper case before
    lower case). The default value for \code{genotype}
    is \code{"freq"}.  The default value for \code{haplotype} is
    \code{"no"}.
  }
  \item{allow.partial.missing}{logical indicating whether one allele is
    permitted to be missing.  When set to \code{FALSE} both alleles
    are set to \code{NA} when either is missing.}
  \item{locus}{ object of class locus, gene, or marker, holding
    information about the source of this genotype.}
  \item{genotypeOrder}{character, vector of genotype/haplotype names so
    that further functions can sort genotypes/haplotypes in wanted
    order}
  \item{...}{optional arguments}
}

\details{

  Genotype objects hold information on which gene or marker alleles were
  observed for different individuals.  For each individual, two alleles
  are recorded.

  The genotype class considers the stored alleles to be unordered, i.e., "C/T"
  is equivalent to "T/C".  The haplotype class considers the order of the
  alleles to be significant so that "C/T" is distinct from "T/C".

  When calling \code{genotype} or \code{haplotype}:

  \itemize{
  \item If only \code{a1} is provided and is a character vector, it is
    assumed that each element encodes both alleles. In this case, if
    \code{sep} is a character string, \code{a1} is assumed to be coded
    as "Allele1<sep>Allele2".  If \code{sep} is a numeric value, it is
    assumed that character locations \code{1:sep} contain allele 1 and
    that remaining locations contain allele 2.

  \item If \code{a1} is a matrix, it is assumed that column 1 contains
  allele 1 and column 2 contains allele 2.

  \item If \code{a1} and \code{a2} are both provided, each is assumed to
  contain one allele value so that the genotype for an individual is
  obtained by \code{paste(a1,a2,sep="/")}.

  }

  If \code{remove.spaces} is TRUE, (the default) any whitespace
  contained in \code{a1} and \code{a2} is removed when the genotypes are
  created.  If whitespace is used as the separator, (eg "C C", "C T",
  ...), be sure to set remove.spaces to FALSE.

  When the alleles are explicitly specified using the \code{alleles}
  argument, all potential alleles not present in the list will be
  converted to \code{NA}.

  NOTE: \code{genotype} assumes that the order of the alleles is not important
  (E.G., "A/C" == "C/A").  Use class \code{haplotype} if order is significant.

  If \code{genotypeOrder=NULL} (the default setting), then
  \code{\link{expectedGenotypes}} is used to get standard sorting order.
  Only unique values in \code{genotypeOrder} are used, which in turns
  means that the first occurrence prevails. When \code{genotypeOrder} is
  given some genotype names, but not all that appear in the data, the
  rest (those in the data and possible combinations based on allele
  variants) is automatically added at the end of
  \code{genotypeOrder}. This puts "missing" genotype names at the end of
  sort order. This feature is especially useful when there are a lot of
  allele variants and especially in haplotypes. See examples.

}

\value{

  The genotype class extends "factor" and haplotype extends
  genotype. Both classes have the following attributes:
  \item{levels}{ character vector of possible genotype/haplotype values
    stored coded by \code{paste( allele1, "/", allele2, sep="")}.}
  \item{allele.names}{ character vector of possible alleles. For a SNP,
    these might be c("A","T").   For a variable length dinucleotyde
    repeat this might be c("136","138","140","148"). }
  \item{allele.map}{ matrix encoding how the factor levels correspond to
    alleles.  See the source code to \code{allele.genotype()} for how to
    extract allele values using this matrix.  Better yet, just use
    \code{allele.genotype()}.}
  \item{genotypeOrder}{character, genotype/haplotype names in
    defined order that can used for sorting in various functions. Note
    that this slot stores both ordered and unordered genotypes
    i.e. "A/B" and "B/A".}
}

%\references{ ~put references to the literature/web site here ~ }
\author{Gregory R. Warnes \email{greg@warnes.net} and
  Friedrich Leisch.}
\seealso{
  \code{\link{HWE.test}},
  \code{\link{allele}},
  \code{\link{homozygote}},
  \code{\link{heterozygote}},
  \code{\link{carrier}},
  \code{\link{summary.genotype}},
  \code{\link{allele.count}},
  \code{\link{sort.genotype}},
  \code{\link{genotypeOrder}},
  \code{\link{locus}},
  \code{\link{gene}},
  \code{\link{marker}}, and
  \code{\link{\%in\%}} for default \%in\% method
  }

\examples{
# several examples of genotype data in different formats
example.data   <- c("D/D","D/I","D/D","I/I","D/D",
                    "D/D","D/D","D/D","I/I","")
g1  <- genotype(example.data)
g1

example.data2  <- c("C-C","C-T","C-C","T-T","C-C",
                    "C-C","C-C","C-C","T-T","")
g2  <- genotype(example.data2,sep="-")
g2


example.nosep  <- c("DD", "DI", "DD", "II", "DD",
                    "DD", "DD", "DD", "II", "")
g3  <- genotype(example.nosep,sep="")
g3

example.a1 <- c("D",  "D",  "D",  "I",  "D",  "D",  "D",  "D",  "I",  "")
example.a2 <- c("D",  "I",  "D",  "I",  "D",  "D",  "D",  "D",  "I",  "")
g4  <- genotype(example.a1,example.a2)
g4

example.mat <- cbind(a1=example.a1, a1=example.a2)
g5  <- genotype(example.mat)
g5

example.data5  <- c("D   /   D","D   /   I","D   /   D","I   /   I",
                    "D   /   D","D   /   D","D   /   D","D   /   D",
                    "I   /   I","")
g5  <- genotype(example.data5,rem=TRUE)
g5

# show how genotype and haplotype differ
data1 <- c("C/C", "C/T", "T/C")
data2 <- c("C/C", "T/C", "T/C")

test1  <- genotype( data1 )
test2  <- genotype( data2 )

test3  <-  haplotype( data1 )
test4  <-  haplotype( data2 )

test1==test2
test3==test4

test1=="C/T"
test1=="T/C"

test3=="C/T"
test3=="T/C"

## also
test1 %in% test2
test1 %in% data2
test3 %in% test4

test1 %in% "C/T"
test1 %in% "T/C"

test3 %in% "C/T"
test3 %in% "T/C"

## "Messy" example

m3  <-  c("D D/\t   D D","D\tD/   I",  "D D/   D D","I/   I",
          "D D/   D D","D D/   D D","D D/   D D","D D/   D D",
          "I/   I","/   ","/I")

genotype(m3)
summary(genotype(m3))

m4  <-  c("D D","D I","D D","I I",
          "D D","D D","D D","D D",
          "I I","   ","  I")

genotype(m4,sep=1)
genotype(m4,sep=" ",remove.spaces=FALSE)
summary(genotype(m4,sep=" ",remove.spaces=FALSE))

m5  <-  c("DD","DI","DD","II",
          "DD","DD","DD","DD",
          "II","   "," I")
genotype(m5,sep=1)
haplotype(m5,sep=1,remove.spaces=FALSE)

g5  <- genotype(m5,sep="")
h5  <- haplotype(m5,sep="")

heterozygote(g5)
homozygote(g5)
carrier(g5,"D")

g5[9:10]  <- haplotype(m4,sep=" ",remove=FALSE)[1:2]
g5

g5[9:10]
allele(g5[9:10],1)
allele(g5,1)[9:10]

# drop unused alleles
g5[9:10,drop=TRUE]
h5[9:10,drop=TRUE]

# Convert allele.counts into genotype

x <- c(0,1,2,1,1,2,NA,1,2,1,2,2,2)
g <- as.genotype.allele.count(x, alleles=c("C","T") )
g

# Use of genotypeOrder
example.data   <- c("D/D","D/I","I/D","I/I","D/D",
                    "D/D","D/I","I/D","I/I","")
summary(genotype(example.data))
genotypeOrder(genotype(example.data))

summary(genotype(example.data, genotypeOrder=c("D/D", "I/I", "D/I")))
summary(genotype(example.data, genotypeOrder=c(              "D/I")))
summary(haplotype(example.data, genotypeOrder=c(             "I/D", "D/I")))
example.data <- genotype(example.data)
genotypeOrder(example.data) <- c("D/D", "I/I", "D/I")
genotypeOrder(example.data)

}
\keyword{ misc }
%%\keyword{genetics}
