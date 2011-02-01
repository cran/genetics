% $Id: makeGenotypes.Rd 1346 2011-02-02 02:21:44Z warnes $

\name{makeGenotypes}
\alias{makeGenotypes}
\alias{makeHaplotypes}
\title{Convert columns in a dataframe to genotypes or haplotypes}
\description{ 
  Convert columns in a dataframe to genotypes or haplotypes.
  }
\usage{
makeGenotypes(data, convert, sep = "/", tol = 0.5, ..., method=as.genotype)
makeHaplotypes(data, convert, sep = "/", tol = 0.9, ...)
}
\arguments{
  \item{data}{Dataframe containing columns to be converted}
  \item{convert}{Vector or list of pairs specifying which columns
    contain genotype/haplotype data.   See below for details.}
  \item{sep}{Genotype separator}
  \item{tol}{See below.}
  \item{\dots}{Optional arguments to as.genotype function}
  \item{method}{Function used to perform the conversion.}
}
\details{
  The functions makeGenotypes and makeHaplotypes allow the conversion of
  all of the genetic variables in a dataset to genotypes or haplotypes
  in a single step.

  The parameter \code{convert} may be missing, a vector of
  column names, indexes or true/false indictators, or a list of column
  name or index pairs.

  When the argument \code{convert} is not provided, the function will
  look for columns where at least \code{tol}*100\% of the records
  contain the separator character \code{sep} ('/' by default).  These
  columns will then be assumed to contain both of the genotype/haplotype
  alleles and will be converted in-place to genotype variables.
  
  When the argument \code{convert} is a vector of column names, indexes
  or true/false indictators, the corresponding columns will be assumed
  to contain both of the genotype/haplotype alleles and will be
  converted in-place to genotype variables.

  When the argument \code{convert} is a list containing column name or
  index pairs, the two elements of each pair will be assumed to contain the
  individual alleles of a genotype/haplotype.  The first column
  specified in each pair will be replaced with the new
  genotype/haplotype variable named \code{name1 + sep + name2}.  The
  second column will be removed.

  Note that the \code{method} argument may be used to supply a
  non-standard conversion function, such as
  \code{as.genotype.allele.count}, which converts from [0,1,2] to
  ['A/A','A/B','A/C'] (or the specified allele names).  See the example
  below.
}
\value{
  Dataframe containing converted genotype/haplotype variables. All other
  variables will be unchanged.
}
\author{ Gregory R. Warnes \email{greg@warnes.net }
  }
\seealso{
  \code{\link{genotype}}
}
\examples{
\dontrun{
# common case
data <- read.csv(file="genotype_data.csv")
data <- makeGenotypes(data)
}

# Create a test data set where there are several genotypes in columns
# of the form "A/T".
test1 <- data.frame(Tmt=sample(c("Control","Trt1","Trt2"),20, replace=TRUE),
                G1=sample(c("A/T","T/T","T/A",NA),20, replace=TRUE),
                N1=rnorm(20),
                I1=sample(1:100,20,replace=TRUE),
                G2=paste(sample(c("134","138","140","142","146"),20,
                                replace=TRUE),
                         sample(c("134","138","140","142","146"),20,
                                replace=TRUE),
                         sep=" / "),
                G3=sample(c("A /T","T /T","T /A"),20, replace=TRUE),
                comment=sample(c("Possible Bad Data/Lab Error",""),20,
                               rep=TRUE)
                )
test1

# now automatically convert genotype columns
geno1 <- makeGenotypes(test1)
geno1

# Create a test data set where there are several haplotypes with alleles
# in adjacent columns.
test2 <- data.frame(Tmt=sample(c("Control","Trt1","Trt2"),20, replace=TRUE),
                    G1.1=sample(c("A","T",NA),20, replace=TRUE),
                    G1.2=sample(c("A","T",NA),20, replace=TRUE),
                    N1=rnorm(20),
                    I1=sample(1:100,20,replace=TRUE),
                    G2.1=sample(c("134","138","140","142","146"),20,
                                replace=TRUE),
                    G2.2=sample(c("134","138","140","142","146"),20,
                                replace=TRUE),
                    G3.1=sample(c("A ","T ","T "),20, replace=TRUE),
                    G3.2=sample(c("A ","T ","T "),20, replace=TRUE),
                    comment=sample(c("Possible Bad Data/Lab Error",""),20,
                                   rep=TRUE)
                   ) 
test2

# specifly the locations of the columns to be paired for haplotypes
makeHaplotypes(test2, convert=list(c("G1.1","G1.2"),6:7,8:9))

# Create a test data set where the data is coded as numeric allele
# counts (0-2).
test3 <- data.frame(Tmt=sample(c("Control","Trt1","Trt2"),20, replace=TRUE),
                    G1=sample(c(0:2,NA),20, replace=TRUE),
                    N1=rnorm(20),
                    I1=sample(1:100,20,replace=TRUE),
                    G2=sample(0:2,20, replace=TRUE),
                    comment=sample(c("Possible Bad Data/Lab Error",""),20,
                                   rep=TRUE)
                   ) 
test3

# specifly the locations of the columns, and a non-standard conversion
makeGenotypes(test3, convert=c('G1','G2'), method=as.genotype.allele.count)


}
\keyword{ misc }
