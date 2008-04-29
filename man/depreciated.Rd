% $Id: power.casectl.Rd 1272 2007-07-18 11:57:28Z ggorjan $
\name{Depreciated}
\alias{power.casectrl}
\title{Depreciated functions}
\description{
  These functions are depreciated.
}
\usage{
power.casectrl(...)
}
\arguments{
  \item{\dots}{All arguments are ignored}
}
\details{
  The \code{power.casectl} function contained serious errors and has been
  replaced by \code{\link[GeneticsDesign]{GPC}},
  \code{\link[GeneticsDesign]{GeneticPower.Quantitative.Factor}}, or
  \code{\link[GeneticsDesign]{GeneticPower.Quantitative.Numeric}} in the
  BioConductor GeneticsDesign package.
      
  In specific, the \code{power.casectl} function used an expected
  contingency table to create the test statistic that was
  erroneously based on the underlying null, rather than on the
  marginal totals of the observed table. In addition, the modeling of
  dominant and recessive modes of inheritance had assumed a "perfect" 
  genotype with no disease, whereas in reality a dominant or recessive
  mode of inheritance simply means that two of the genotypes will have
  an identical odds ratio compared to the 3rd genotype (the other
  homozygote). 
}
\keyword{misc}
