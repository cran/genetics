#
# $Id: HWE.test.R,v 1.13 2002/11/08 17:20:05 warnesgr Exp $
#
# $Log: HWE.test.R,v $
# Revision 1.13  2002/11/08 17:20:05  warnesgr
# - Added test for more than two alleles to diseq.ci.  If there are more
#   than two, generate a warning message indicating that confidence
#   intervals for values near 0 and 1 are ill-behaved.  When
#   'correct=TRUE' use ci.balance to attempt to correct for this.
# - Modified HWE.test to display warnings from diseq.ci.
#
# Revision 1.12  2002/10/28 18:19:03  warnesgr
#
# - Lots of changes to the disequlibrim computation and confidence interval code.
#    - D, D', and R are now computed
#    - Overall values are reported *with sign* of the number of alleles is 2
#    - Bootstrapping code is now more efficient
#
# - Changes to HWE.test() and print.hwe.test() corresponding to the above changes
#
# Revision 1.11  2002/10/11 18:01:00  warnesgr
#
# Updated version number.
#
# Revision 1.10  2002/09/24 01:32:19  warnesgr
# - 'Un-genericized' diseq()
# - Moved documentation of diseq() and diseq from HWE.test.Rd to diseq.Rd
# - Cleaned up HWE.test.Rd and diseq.Rd
# - Fixed broken class check in diseq() and diseq.ci()
# - Removed allele.count.default() -- this will force the user to
#   explicitly call 'genotype' on the data to use allele.count().
# - Added zzz.R to require package 'boot'
#
# Revision 1.9  2002/09/24 00:02:01  WarnesGR
# - Moved code that computed D-hat to diseq() and diseq.genotype() in diseq.R.
# - Added diseq.ci() to compute confidence interval for D-hat.
# - Added code to call diseq() and diseq.ci() from HWE.test()
# - Added arguments to HWE.test() and print.HWE.test() to control these new features
# - Added text to HWE.test.Rd documenting these new functions and arguments
#
# Revision 1.8  2002/06/27 18:46:04  warnesgr
# - Allow user to specify parameters for the chisquare test.
#
# Revision 1.7  2002/06/19 00:13:54  warnesgr
# - Fixed bug that caused 'Error: subscript out of bounds' when a 2
#   allele genotype had 1/1 and 1/2 but not 2/2.
#
# Revision 1.6  2002/02/14 12:55:14  warnes
#
# Fixed bugs in [.genotype and [.haplotype.
#
# - There was a serious typo that would have been cought had I just done
#   R CMD check.
#
# - Revamped [.genotype and [.haplotype to work correctly and to avoid
#   the overhead of calling genotype() or haplotype() when drop=TRUE.  I
#   think these finally work properly.
#
# Revision 1.5  2001/08/13 15:36:45  warnes
#
# BUGFIX: Corrected problem where the allele order was not consistent between
# 	the first (Observed #) and second columns (Expected #) of the
# 	HWE.test table.  This caused incorrect values for Obs-Exp to
# 	be calcuated, as well as incorrect D-hat estimates.
# 	Fortunately, it did not effect the correctness of the test
# 	statistic or p-value, since these were computed via the chisq
# 	function.
#
# Revision 1.4  2001/06/28 19:14:00  warnes
#
# - Modified to use chisq.test(...,simulate.p.value=TRUE,B=10000) to use
# simulated p-value so that assumptions underlying Chi-square
# approximation need not hold.  This particularly useful for small
# sample sizes or when some allele pairs are not observed.
#
# - Now gives a table of observed counts, observed frequency, expected
# frequency, and disequilibrium statistic (D-hat) for each pair of
# alleles in addition to the overall disequilibrium statistic (D-hat).
#
# Revision 1.3  2001/06/15 17:02:19  warnes
#
# - Modified HWE.test to use chisq.test.  Added a modified version of
#   chisq.test until the (minor) changes show up in ctest.
#
# - Fixed documentation links and alias omissions
#
# - Fixed getallele() to snag allele element or attribute if either is
#   present.  This makes it work for summary.genotype as well as for
#   genotype.
#
# Revision 1.2  2001/06/13 17:03:33  warnes
#
# - Updated HWE.test.R to not assume that a summary.genotype object
# contains the original genotype object.
#
# - Removed a duplicate link to "allele.names" from genotype.Rd.
#
# Revision 1.1  2001/05/30 22:12:34  warnes
# Updated documentation mostly.  Added as.character.genotype().
#
#


####
### Hardy-Weinburg Equilibrium Test
###

HWE.test <- function(x, ...)
{
	UseMethod("HWE.test")
}


HWE.test.genotype <- function(x, simulate.p.value=TRUE, B=10000,
                              conf=0.95, ci.B=1000, ... )
  # future options "bootstrap","exact"
{

  retval <- list()

  # compute disequlibrium
  retval$diseq <- diseq(x)
  
  # compute confidence intervals
  retval$ci <- diseq.ci(x, R=ci.B, conf=conf)


  # compute p-value
  tab <- retval$diseq$observed.no
  tab  <- 0.5 * (tab + t(tab))   # make symmetric for chisq.test

  retval$test  <- chisq.test(tab,simulate.p.value=simulate.p.value,B=B,...
                             )
                             # , df= k*(k-1)/2 ) # when possible


  retval$simulate.p.value <- simulate.p.value
  retval$B <- B
  retval$conf <- conf
  retval$ci.B <- ci.B
  retval$test$data.name  <- deparse(substitute(x))
  retval$call  <- match.call()
  class(retval)  <-  c("HWE.test")
  return(retval)
}



print.HWE.test  <-  function(x, show=c("D","D'","r"), ...)
  {

    cat("\n")
    cat("\t-----------------------------------\n")
    cat("\tTest for Hardy-Wienburg-Equilibrium\n")
    cat("\t-----------------------------------\n")
    cat("\n")
    if(!is.null(x$locus))
      {
        cat("\n")
        print( x$locus )
      }
    cat("Call: \n")
    print(x$call)
    cat("\n")
    if("D" %in% show)
      {
        cat("Raw Disequlibrium for each allele pair (D)\n")
        cat("\n") 
        print(x$diseq$D)
        cat("\n")
      }
    if("D'" %in% show)
      {
        cat("Scaled Disequlibrium for each allele pair (D')\n")
        cat("\n") 
        print(x$diseq$Dprime)
        cat("\n")
      }
    if("r" %in% show)
      {
        cat("Correlation coefficient for each allele pair (r)\n")
        cat("\n") 
        print(x$diseq$r)
        cat("\n")
      }
    cat("Overall Values (mean absolute-value weighted by expected allele frequency)\n")
    cat("\n")

    show.tab <- NULL
    
    if("D" %in% show)
      show.tab <- rbind(show.tab, "  D"=x$diseq$D.overall)
    if("D'" %in% show)
      show.tab <- rbind(show.tab, "  D'"=x$diseq$Dprime.overall)
    if("r" %in% show)
      show.tab <- rbind(show.tab, "  r"=x$diseq$r.overall)

    colnames(show.tab) <- "Value"

    print(show.tab)
    
    cat("\n") 

    whichvec <- c("D","D'","r") %in% show

    cat("Confidence intervals computed via bootstrap using", x$ci.B, "samples\n")
    cat("\n")

    if(!is.null(x$ci$warning.text))
      cat(strwrap(paste("WARNING:", x$ci$warning.text), prefix="    * "),"\n",
          sep="\n")
    
    show.tab <- matrix(ncol=4, nrow=3)
    tmp <- format(x$ci$ci[,1:3], digits=getOption("digits"))
    show.tab[,1] <- tmp[,1]  # Observed
    show.tab[,2] <- paste("(", tmp[,2], ", ", tmp[,3], ")", sep="" )
    show.tab[,3] <- x$ci$ci[,4]
    show.tab[,4] <- ifelse(x$ci$ci[,5],"YES","*NO*")

    colnames(show.tab) <- c("Observed", "95% CI", "NA's", "Contains Zero?")
    rownames(show.tab) <- paste("  ", rownames(tmp), sep="")
    
    print(show.tab[whichvec,], quote=FALSE)

    cat("\n")
    cat("Significance Test:\n")
    print(x$test)
    cat("\n")
    cat("\n")
  }

