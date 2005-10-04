# expectedGenotypes.R
#--------------------------------------------------------------------------
# What: Construct expected genotypes according to known allele variants
# Time-stamp: <2005-03-08 11:43:18 ggorjan>
#--------------------------------------------------------------------------

expectedGenotypes <- function (x, alleles=allele.names(x), ploidy=2) {
    require("gtools")
    # Checks
    if (missing(x) && missing(alleles))  
        stop("at least one of 'x' or 'alleles' must be given")
    if (!(missing(x) && !missing(alleles))) {
        if (!is.genotype(x)) stop("x must be of class 'genotype'")
    }
    # Find possible genotypes according to allele variants
    comb <- combinations(n=length(alleles),
                         r=ploidy,
                         v=alleles,
                         repeats.allowed=T)
    # Create a nice character vector of expected genotypes
    genotypes <- vector(mode="character", length=dim(comb)[1])
    for (i in 1:dim(comb)[1]) {
        tmp1 <- as.character(comb[i,1])
        if (ploidy >= 2) {
            for (j in 2:ploidy) {
                tmp1 <- paste(tmp1, "/", comb[i,j], sep="")
            }
        }
        genotypes[i] <- tmp1
    }
    return(genotypes)
}

#--------------------------------------------------------------------------
# expectedGenotypes.R ends here
