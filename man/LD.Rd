% $Id: LD.Rd 1346 2011-02-02 02:21:44Z warnes $

\name{LD}
\alias{LD}
\alias{LD.genotype}
\alias{LD.data.frame}
\title{Pairwise linkage disequilibrium between genetic markers.}
\description{
  Compute pairwise linkage disequilibrium between genetic markers
}
\usage{
LD(g1, ...)
\method{LD}{genotype}(g1,g2,...)
\method{LD}{data.frame}(g1,...)
}
\arguments{
  \item{g1}{ genotype object or dataframe containing genotype objects }
  \item{g2}{ genotype object (ignored if g1 is a dataframe) }
  \item{\dots}{ optional arguments (ignored) }
}
\details{
  Linkage disequilibrium (LD) is the non-random association of
  marker alleles and can arise from marker proximity or from selection
  bias.
  
  \code{LD.genotype} estimates the extent of LD for a single pair of
  genotypes.  \code{LD.data.frame} computes LD for all pairs of
  genotypes contained in a data frame.  Before starting,
  \code{LD.data.frame} checks the class and number of alleles of each
  variable in the dataframe.  If the data frame contains non-genotype
  objects or genotypes with more or less than 2 alleles, these will be
  omitted from the computation and a warning will be generated.

  Three estimators of LD are computed:
  \itemize{

    \item{D}{ raw difference in frequency between the
              observed number of AB pairs and the expected number:
      
      \deqn{%
 	D = p_{AB} - p_A p_B %
      }{%
 	D = p(AB) - p(A)*p(B) %
      }
      
    }
    \item{D'}{ scaled D spanning the range [-1,1] 
      
      
      \deqn{D' = \frac{D}{D_{max} } }{D' = D / Dmax}
      
      where, if D > 0:
      \deqn{%
  	D_{max} = \min( p_A p_b, p_a p_B )  %
      }{%
  	Dmax = min( p(A)p(b), p(a)p(B) )   %
      } 
      or if D < 0:
      \deqn{%
  	D_{max} = \max{ -p_A p_B, -p_a p_b }  %
      }{%
  	Dmax = max( -p(A)p(B), -p(a)p(b) )  %
      }
    }
    
    \item{r}{ correlation coefficient between the markers
      
      \deqn{%
 	r = \frac{-D}{\sqrt( p_A * p_a * p_B * p_b  )} %
      }{%
 	r = -D / sqrt( p(A) * p(a) * p(B) * p(b) ) %
      }
    }
  }
  
  where
  \itemize{
    \item{-}{ \eqn{p_A}{p(A)} is defined as the observed probability of
      allele 'A' for marker 1, }
    \item{-}{ \eqn{p_a=1-p_A}{p(a) = 1-p(A)} is defined as the observed probability of
      allele 'a' for marker 1, }
    \item{-}{\eqn{p_B}{p(B)} is defined as the observed probability of
      allele 'B' for marker 2, and }
    \item{-}{\eqn{p_b=1-p_B}{p(b) = 1- p(B)} is defined as the observed probability of
      allele 'b' for marker 2, and }
    \item{-}{\eqn{p_{AB}}{p(AB)} is defined as the probability of
      the marker allele pair 'AB'. }
  }

  For genotype data, AB/ab cannot be distinguished from
  aB/Ab. Consequently, we estimate \eqn{p_{AB}}{p(AB)} using maximum
  likelihood and use this value in the computations.
  }
\value{

  \code{LD.genotype} returns a 5 element list:
    \item{call}{the matched call}
    \item{D}{Linkage disequilibrium estimate}
    \item{Dprime }{Scaled linkage disequilibrium estimate}
    \item{corr}{Correlation coefficient} 
    \item{nobs}{Number of observations}
    \item{chisq}{Chi-square statistic for linkage
      equilibrium (i.e., D=D'=corr=0)}
    \item{p.value}{Chi-square p-value for marker independence}

  \code{LD.data.frame} returns a list with the same elements, but each
  element is a matrix where the upper off-diagonal elements contain the
  estimate for the corresponding pair of markers.  The other matrix
  elements are \code{NA}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Gregory R. Warnes \email{greg@warnes.net} }
\seealso{ \code{\link{genotype}}, \code{\link{HWE.test}}  }
\examples{

g1 <- genotype( c('T/A',    NA, 'T/T',    NA, 'T/A',    NA, 'T/T', 'T/A',
                  'T/T', 'T/T', 'T/A', 'A/A', 'T/T', 'T/A', 'T/A', 'T/T',
                     NA, 'T/A', 'T/A',   NA) )

g2 <- genotype( c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
                  'C/A', 'C/C', 'C/A', 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
                  'C/A', 'C/A', 'C/A', 'A/A') )


g3 <- genotype( c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
                  'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A', 'T/T',
                  'T/A', 'T/A', 'T/A', 'T/T') )

# Compute LD on a single pair

LD(g1,g2)

# Compute LD table for all 3 genotypes

data <- makeGenotypes(data.frame(g1,g2,g3))
LD(data)
}
\keyword{misc}
