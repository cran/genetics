# $Id: HWE.exact.R,v 1.2 2002/12/02 16:12:00 warnesgr Exp $
#
# $Log: HWE.exact.R,v $
# Revision 1.2  2002/12/02 16:12:00  warnesgr
#
# - Now returnes only the computed P-value.
#
# Revision 1.1  2002/12/02 15:40:40  warnesgr
#
# - Added code for HWE.exact, based on contributed function by David
#   Duffy <davidD@qimr.edu.au>
#
#
# Based on code submitted by David Duffy <davidD@qimr.edu.au>
#
# Exact test for HWE: 2 alleles
#
# See eg Emigh TH.  Comparison of tests for Hardy-Weinberg Equilibrium.
#  Biometrics 1980; 36: 627-642
#

HWE.exact <- function(x)
{
  if(!is.genotype(x))
    stop("x must be of class 'genotype' or 'haplotype'")

  nallele <- length(na.omit(allele.names(x)))

  if(nallele != 2)
    stop("Exact HWE test can only be computed for 2 markers with 2 alleles")


  allele.tab <- table( factor(allele(x,1), levels=allele.names(x)),
                       factor(allele(x,2), levels=allele.names(x)) )

  n11 <- allele.tab[1,1]
  n12 <- allele.tab[1,2] + allele.tab[2,1]
  n22 <- allele.tab[2,2]
  
  n1 <- 2*n11+n12
  n2 <- 2*n22+n12

  
  dhwe2 <- function(n11, n12, n22) {
    f <- function(x) lgamma(x+1)
    n <- n11+n12+n22
    n1 <- 2*n11+n12
    n2 <- 2*n22+n12
    exp(log(2)*(n12) + f(n) - f(n11) - f(n12) - f(n22) - f(2*n) +
        f(n1) + f(n2))
  }


  x12 <- seq(n1 %% 2,min(n1,n2),2)
  x11 <- (n1-x12)/2
  x22 <- (n2-x12)/2
  dist <- data.frame(n11=x11,n12=x12,n22=x22,density=dhwe2(x11,x12,x22))
  dist <- dist[order(dist$density),]
  return(prob=cumsum(dist$density)[dist$n11==n11 &
                                       dist$n12==n12 &
                                       dist$n22==n22])
}

