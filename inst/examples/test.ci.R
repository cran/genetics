# test 3 allele model

library(genetics)
library(combinat)

gen3 <- function(d=0, nobs=20)
  {
    pvec <- c(aa=1+d, ab=2*(1-d), ac=2*(1-d), bb=1+d, bc=2*(1-d), cc=1+d)
    pvec <- pvec/sum(pvec)

    gen <- rmultinomial(n=nobs, p=rbind(pvec), rows=1)
    genotype( rep( c("A/A","A/B","A/C","B/B","B/C","C/C"), gen), alleles=c("A","B","C") )
  }

worker <- function(...)
  diseq.ci(gen3())$ci

sim <- t(sapply(1:10, function(x) worker))
