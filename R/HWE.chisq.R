# $Id: HWE.chisq.R,v 1.2 2003/05/22 17:25:23 warnesgr Exp $


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

