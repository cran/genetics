#
# $Id: HWE.chisq.R,v 1.1 2003/03/07 14:32:48 warnesgr Exp $
#
# $Log: HWE.chisq.R,v $
# Revision 1.1  2003/03/07 14:32:48  warnesgr
#
# - Created HWE.chisq, HWE.chisq.genotype and corresponding man page.
# - Moved computation of Chisquare test for HWE from HWE.test.genotype
#   to HWE.chisq.genotype.
# - Added option (on by default) to compute the exact p-value using HWE.exact.
#   This is on by default when nallele=2
#
#

###
### Hardy-Weinberg Equilibrium Significance Test
###


HWE.chisq <- function(x, ...)
  UseMethod("HWE.chisq")

HWE.chisq.genotype <- function (x, simulate.p.value = TRUE, B = 10000, ...)
{
    observed.no <- table(factor(allele(x, 1), levels = allele.names(x)), 
        factor(allele(x, 2), levels = allele.names(x)))
    tab <- observed.no
    tab <- 0.5 * (tab + t(tab))
    k <- ncol(tab)
    if(simulate.p.value)
      {
        test <- chisq.test(tab, simulate.p.value = simulate.p.value, 
                           B = B, ...)
      }
    else
      {
        test <- chisq.test(tab, ...)
        test$parameter <- k*(k-1)/2
        test$p.value <- pchisq(test$statistic, test$parameter, lower = FALSE)
        names(test$statistic) <- "X-squared"
        names(test$parameter) <- "df"
      }
    return(test)
}

