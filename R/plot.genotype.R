## plot.genotype.R
###------------------------------------------------------------------------
## What: Plot genotype object
## $Id: plot.genotype.R,v 1.4 2005/10/04 16:30:27 warnes Exp $
## Time-stamp: <2005-09-13 22:11:19 ggorjan>
###------------------------------------------------------------------------

plot.genotype <- function(x,
                          type=c("genotype", "allele"),
                          what=c("percentage","number"),
                          ...)
{
  what <- match.arg(what)
  type <- match.arg(type)

  ## get details
  tmp <- summary(x)
  
  ## Percentages or numbers
  if (what == "percentage")
    whati <- 2
  else
    whati <- 1
  
  ## Plot
  if (type == "allele")
    barplot(tmp$allele.freq[, whati], ...)
  else # genotype
    barplot(tmp$genotype.freq[, whati], ...)
}

###------------------------------------------------------------------------
## plot.genotype.R ends here
