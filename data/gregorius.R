# $Id: gregorius.R,v 1.1 2002/12/02 17:13:29 warnesgr Exp $
#
# $Log: gregorius.R,v $
# Revision 1.1  2002/12/02 17:13:29  warnesgr
#
# - Added data set (and documentation) 'gregorius' which contains
#   necessary samples sizes to detect alleles.
#
#
#                                        
# Power to detect the range of alleles in a population
#
# Sample sizes required to identify *all* alleles with frequency >= alpha
# with a probability >=sigma
#
# From:
#
# Gregorius, H.-R. 1980. The probability of losing an allele when
# diploid genotypes are sampled.  Biometrics 36, 643-652.
#
# Submitted by David Duffy <davidD@qimr.edu.au>
#

gregorius <- data.frame(
  minfreq=c(0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05, 0.04,0.03,
            0.02,0.01, 0.009, 0.008),
  sigma05=c(6,7,11,21,51,57,65,77,92,117,152,212,341,754,850,972),
  sigma01=c(8,10,15,28,66,74,84,99,119,149,192,265,422,916,1030,1174),
  sigma001=c(11,14,22,39,88,99,112,131,156,194,249,341,536,1146,1285,1462)
)

