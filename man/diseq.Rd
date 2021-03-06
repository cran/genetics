% $Id: diseq.Rd 1346 2011-02-02 02:21:44Z warnes $

\name{diseq}
\alias{diseq}
\alias{diseq.table}
\alias{diseq.genotype}
\alias{diseq.ci}
\alias{print.diseq}
\title{Estimate or Compute Confidence Interval for the Single-Marker Disequilibrium}
\description{
  Estimate or compute confidence interval for single-marker disequilibrium.
}
\usage{
diseq(x, ...)
\method{print}{diseq}(x, show=c("D","D'","r","R^2","table"), ...)
diseq.ci(x, R=1000, conf=0.95, correct=TRUE, na.rm=TRUE, ...)
}
\arguments{
  \item{x}{genotype or haplotype object.}
  \item{show}{a character value or vector indicating which
    disequilibrium measures should be displayed.  The default is to show
    all of the available measures. \code{show="table"} will display
    a table of observed, expected, and observed-expected
    frequencies.}
  \item{conf}{Confidence level to use when computing the confidence
    level for D-hat.  Defaults to 0.95, should be in (0,1). }
  \item{R}{Number of bootstrap iterations to use when computing the
    confidence interval. Defaults to 1000.}
  \item{correct}{See details.}
  \item{na.rm}{ logical. Should missing values be removed?}
  \item{...}{ optional parameters passed to \code{boot.ci}
    (\code{diseq.ci}) or ignored.}
}
\details{
  For a single-gene marker, \code{diseq} computes the Hardy-Weinberg
  (dis)equilibrium statistic D, D', r (the correlation coefficient), and
  \eqn{r^2} for each pair of allele values, as well as an overall
  summary value for each measure across all alleles.  \code{print.diseq}
  displays the contents of a \code{diseq} object. \code{diseq.ci}
  computes a bootstrap confidence interval for this estimate.

  For consistency, I have applied the standard definitions for D, D',
  and r from the Linkage Disequilibrium case, replacing all marker 
  probabilities with the appropriate allele probabilities.
  
  Thus, for each allele pair,
  \itemize{
    \item{D}{ is defined as the half of the raw difference
      in frequency between
      the observed number of heterozygotes and the expected number:
      
      \deqn{%
 	D = \frac{1}{2} ( p_{ij} + p_{ji} ) - p_i p_j  %
      }{%
 	D = 1/2 * ( p(ij) + p(ji) ) - p(i)*p(j)   %
      }
      
    }
    \item{D'}{ rescales D to span the range [-1,1] 
      
      
      \deqn{D' = \frac{D}{D_{max} } }{D' = D / Dmax}
      
      where, if D > 0:
      \deqn{%
  	D_{max} = \min{ p_i p_j, p_j p_i } = p_i p_j   %
      }{%
  	Dmax = min(p(i)p(j), p(j)p(i)) =  p(i)p(j)   %
      } 
      or if D < 0:
      \deqn{%
  	D_{max} = \min{ p_i (1 - p_j), p_j (1 - p_i) } %
      }{%
  	Dmax = min( p(i) * (1 - p(j)), p(j)( 1 - (1-p(i) ) ) )
      }
    }
    
    \item{r}{ is the correlation coefficient between two alleles,
      %ignoring all other alleles,
      and can be computed by
      
      \deqn{%
 	r = \frac{-D}{\sqrt( p_i * (1-p_i) p(j) (1-p_j ) )} %
      }{%
 	r = -D / sqrt( p(i)*(1-p(i)) * p(j)*(1-p(j)) ) %
      }
    }
    
  }
  
  where
  \itemize{
    \item{-}{ \eqn{p_i}{p(i)} defined as the observed probability of
      allele 'i', }
    \item{-}{\eqn{p_j}{p(j)} defined as the observed probability of
      allele 'j', and }
    \item{-}{\eqn{p_{ij}}{p(ij)} defined as the observed probability of
      the allele pair 'ij'. }
  }
  
  When there are more than two alleles, the summary values for these
  statistics are obtained by computing a weighted average of the
  absolute value of each allele pair, where the weight is determined by
  the expected frequency. For example:

   \deqn{%
     D_{overall} = \sum_{i \ne j}  |D_{ij}| * p_{ij} % 
   }{%
     D.overall = sum |D(ij)| * p(ij) %
   }

  Bootstrapping is used to generate confidence interval in order to
  avoid reliance on parametric assumptions, which will not hold for
  alleles with low frequencies (e.g. \eqn{D'} following a a Chi-square 
  distribution).  

  See the function \code{\link[genetics]{HWE.test}} for testing
  Hardy-Weinberg Equilibrium, \eqn{D=0}.
  
}
\value{
  \code{diseq} returns an object of class \code{diseq} with components
  \itemize{
    \item{call}{function call used to create this object}
    \item{data}{2-way table of allele pair counts}
    \item{D.hat}{matrix giving the observed count, expected count,
      observed - expected difference, and estimate of disequilibrium for
      each pair of alleles as well as an overall disequilibrium value.}
    \item{TODO}{more slots to be documented}
  }

  \code{diseq.ci} returns an object of class \code{\link[boot]{boot.ci}}
}
\author{ Gregory R. Warnes \email{greg@warnes.net }
  }

\seealso{
  \code{\link{genotype}},
  \code{\link{HWE.test}},
  \code{\link[boot]{boot}},
  \code{\link[boot]{boot.ci}}
}

\examples{
\testonly{
set.seed(7981357)
}
example.data   <- c("D/D","D/I","D/D","I/I","D/D",
                    "D/D","D/D","D/D","I/I","")
g1  <- genotype(example.data)
g1

diseq(g1)
diseq.ci(g1)
HWE.test(g1)  # does the same, plus tests D-hat=0

three.data   <- c(rep("A/A",8),
                  rep("C/A",20),
                  rep("C/T",20),
                  rep("C/C",10),
                  rep("T/T",3))

g3  <- genotype(three.data)
g3

diseq(g3)
diseq.ci(g3, ci.B=10000, ci.type="bca")

# only show observed vs expected table
print(diseq(g3),show='table')

}
\keyword{misc}
%%\keyword{genetics}
