\name{homozygote}
\alias{homozygote}
\alias{heterozygote}
\alias{carrier}
\alias{allele}
\alias{allele.count}
\alias{allele.names}
\title{Extract Features of Genotype objects}
\description{
  \code{homozygote} creates an vector of logicals that are true when the
  alleles of the corresponding observation are the identical.
  
  \code{heterozygote } creates an vector of logicals that are true when the
  alleles of the corresponding observation differ.
  
  \code{carrier} create a logical vector or matrix of logicals
  indicating whether the specified alleles are present.
  
  \code{allele.count} returns the number of copies of the specified
  alleles carried by each observation.
  
  \code{allele} extract the specified allele(s) as a character vector
  or a 2 column matrix.

  \code{allele.names} extract the set of allele names.
 }
\usage{
homozygote(x, allele.name,...)
heterozygote(x,allele.name,...)
\method{carrier}{genotype}(x, allele.name=allele.names(x), any=!missing(allele),
        na.rm=FALSE, ...)
allele.count(x, allele.name=allele.names(x),any=!missing(allele),
             na.rm=FALSE)
allele(x, which=c(1,2) )
allele.names(x)
}
\arguments{
  \item{x}{ \code{genotype} object }
  \item{\dots}{ optional parameters (ignored) }
  \item{allele.name}{ character value or vector of allele names}
  \item{any}{ logical value.  When \code{TRUE}, a single count or
    indicator is returned by combining the results for all of the
    elements of \code{allele}. If \code{FALSE} separate counts or
    indicators should be returned for each element of
    \code{allele}.  Defaults to \code{FALSE} if \code{allele} is
    missing. Otherwise defaults to \code{TRUE}.}
  \item{na.rm}{ logical value indicating whether to remove missing
    values.  When true, any \code{NA} values will be replaced by
    \code{0} or \code{FALSE} as appropriate.  Defaults to \code{FALSE}.}
  \item{which}{ selects which allele to return. For first allele use
    \code{1}.  For second allele use \code{2}.  For both (the default)
    use \code{c(1,2)}.}
      
}

\details{
  When the \code{allele.name} argument is given, heterozygote and
  homozygote return \code{TRUE} if \emph{exactly} one or both alleles,
  respectively, match the specified allele.name.
}

\value{
  \code{homozygote} and \code{heterozygote } return a vector of
  logicals.

  \code{carrier} returns a logical vector if only one allele is
  specified, or if \code{any} is \code{TRUE}.  Otherwise, it returns
  matrix of logicals with one row for each element of \code{allele}.

  
  \code{allele.count} returns a  vector of counts if only one allele is
  specified, or if \code{any} is \code{TRUE}.  Otherwise, it returns
  matrix of counts with one row for each element of \code{allele}.
  
  \code{allele} returns a character vector when one allele is
  specified.  When 2 alleles are specified, it returns a 2 column
  character matrix.

  \code{allele.names} returns a character vector containing the set of
  allele names.
  
  
}
\author{ Gregory R. Warnes \email{Gregory\_R\_Warnes@groton.pfizer.com} }

\seealso{
  \code{\link{genotype}},
  \code{\link{HWE.test}},
%  \code{\link{allele}},
%  \code{\link{homozygote}},
%  \code{\link{heterozygote}}, 
%  \code{\link{carrier}},
  \code{\link{summary.genotype}},
  %  \code{\link{allele.count}}
  \code{\link{locus}}
  \code{\link{gene}}
  \code{\link{marker}}
  }


  \examples{

example.data   <- c("D/D","D/I","D/D","I/I","D/D","D/D","D/D","D/D","I/I","")
g1  <- genotype(example.data)
g1

heterozygote(g1)
homozygote(g1)

carrier(g1,"D")
carrier(g1,"D",na.rm=TRUE)

# get count of one allele 
allele.count(g1,"D")

# get count of each allele
allele.count(g1)  # equivalent to
allele.count(g1, c("D","I"), any=FALSE)

# get combined count for both alleles
allele.count(g1,c("I","D"))

# get second allele
allele(g1,2)

# get both alleles
allele(g1)

}
\keyword{ misc }
%\keyword{genetics}%-- one or more ...
