# plot.genotype.R
#--------------------------------------------------------------------------
# What: Plot genotype object
# Time-stamp: <2005-03-08 11:43:18 ggorjan>
#--------------------------------------------------------------------------

plot.genotype <- function(x, type=c("genotype","allele"),
                          what=c("percentage","number"),  ...)
{

  type <- match.arg(type)
  what <- match.arg(what)
  
  tmp <- summary(x)
  # Percentages or numbers
  if (what == "percentage")
    whati <- 2
  else 
    whati <- 1
  # Plot
  if (type == "allele")
      barplot(tmp$allele.freq[, whati], ...)
  else # genotype
      barplot(tmp$genotype.freq[, whati], ...)    
}

#--------------------------------------------------------------------------
# plot.genotype.R ends here
