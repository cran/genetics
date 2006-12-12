# $Id: HWE_Test.R 61 2002-06-27 18:46:05Z warnesgr $
#
# $Log$
# Revision 1.2  2002/06/27 18:46:05  warnesgr
# - Allow user to specify parameters for the chisquare test.
#
# Revision 1.1  2001/05/07 13:22:39  warnes
#
# Added example files, code, and output.
#
# Revision 1.4  2001/05/01 14:33:19  warneg
#
# Updated files to use changed PG database output format.  The new format is
#
#    Patient ID,Gene,Marker,Allele1/Allele2
#
# Before it was
#
#    Patient ID,Gene,Marker,Count of Allele1,Count of Allele2,Count of Allele 3, ...
#
# This involved changes in  Allele_Freq.R, HWE_Test.R, and test.data.txt
#
# ---
#
# Modified Examples.R to remove random values.  This will allow
# diffing current and previous versions of the code to check for
# regressions.
#
# ---
#
# Fixed as.genotype.allele.count() to handle conversions both when when NA values
# are and are not obtained.
#
# Simplified the class type of HWE.test results to "HWE.test" from
# "HWE.test.allele.freq".
#
# Revision 1.3  2001/04/25 17:45:37  warneg
# Fixed typo that caused an error.
#
# Revision 1.2  2001/04/23 19:39:01  warneg
# Updated to use revised Genomics.R that provides "genotype" and "haplotype" classes.
#
# Revision 1.1  2001/02/06 23:09:44  warneg
#
#
# HWE_Test.R performs the Hardy-Weinberg equilibrium test for the markers
# supplied in the input file.  Initial revision.
#
# Revision 1.2  2001/02/06 17:00:26  warneg
#
#
# Added CVS tags to track version.
#
#
# run as
#   /usr/local/bin/R --vanilla --slave < Allele<-Freq.R

# first, get the library functions
library(genetics)


# get the name of the file containing the allele data
file.name <- Sys.getenv("ALLELE_INPUT_FILENAME")
if(file.name=="")
  {
  warning(paste("Unable to read input file name from the environment\n",
                "variable '\$ALLELE_INPUT_FILENAME'. ", 
                "Using 'input.data.txt' instead.\n",sep=""));
  file.name <- "test.data.txt"
  }

# get the data
cat("\nReading data file '", file.name, "' ...", "\n", sep="" )
input.data <- read.table(file.name,sep=", ", header=T)

# report on what we have
cat( dim(input.data)[1], " rows and ", dim(input.data)[2], " columns were read. \n\n")

cat("Column names are: ", names(input.data), "\n" )
cat("Note: Spaces and '<-' characters are converted to periods ('.') \n")

# make all names uppercase
names(input.data) <- toupper(names(input.data)) 

# check that we have "PATIENT.ID", "LOCUS", and "MARKER" fields.
# If not give warning and assume these are columns 1, 2, and 3.
if( is.na(match("PATIENT.ID", names(input.data) ) ) )
{
	warning(paste( "No column labeled 'PATIENT ID'.\n",
                       "Assuming that the first column ('",
                       names(input.data)[1],
                       "' contains patient id. \n", sep='') ) 
		
	names(input.data)[1] <- "PATIENT.ID"
}



if( is.na(match("LOCUS", names(input.data) ) ) )
{
	warning(paste( "No column labeled 'LOCUS'.\n",
                       "Assuming that the second column ('",
                       names(input.data)[2],
                       "' contains locus/gene name. \n", sep='')) 
		
	names(input.data)[2] <- "LOCUS"
}
if( is.na(match("MARKER", names(input.data) ) ) )
{
	warning(paste( "No column labeled 'MARKER'.",
                       "Assuming that the third column ('",
                       names(input.data)[3],
                       "' contains marker name. \n", sep='')) 
		
	names(input.data)[3] <- "MARKER"
}

if( is.na(match("GENOTYPE", names(input.data) ) ) )
{
	warning(paste( "No column labeled 'GENOTYPE'.",
                       "Assuming that the fourth column ('",
                       names(input.data)[4],
                       "' contains genotype. \n", sep='')) 
		
	names(input.data)[4] <- "GENOTYPE"
}

#
# convert data to 1 record per patient
#
input.data$LOCUS.MARKER  <-  paste(input.data$LOCUS,input.data$MARKER,sep=":")
data  <- data.frame(PATIENT.ID=unique(as.character(input.data$PATIENT.ID)))
data[,unique(input.data$LOCUS.MARKER)]  <- NA
data  <- as.matrix(data)
rownames(data)  <-  data[,1]

tmp  <- split(input.data[,c("PATIENT.ID","LOCUS.MARKER","GENOTYPE")], input.data$LOCUS.MARKER)

for(i in 1:nrow(input.data))
  data[ as.character(input.data[i,"PATIENT.ID"]), 
       as.character(input.data[i,"LOCUS.MARKER"]) ]  <-
  as.character(input.data[i,"GENOTYPE"])

data  <- data.frame(apply( data[,-1], 2, as.character ))
data  <- data.frame(sapply( data, as.genotype, simplify=F))

## Now iterate through doing the HWE test and displaying output

ind  <- !duplicated(input.data$LOCUS.MARKER)
namemat  <- input.data[ind,c("LOCUS","MARKER","LOCUS.MARKER")]
nmarker <-  sum(ind)

for(i in 1:nmarker)
{
  gene  <- as.character(namemat[i,"LOCUS"])
  marker  <- as.character(namemat[i,"MARKER"])

  cat("\n")
  cat("+-------------------------------------\n");
  if(!is.null(gene))
    cat("|\tGene:\t ", gene, "\n");
  
  if(!is.null(marker))
    cat("|\tMarker:\t ", marker, "\n");
  cat("+-------------------------------------\n");

  # compute and print the allele and genotype frequencies
  sum  <-  summary(data[,i])
  print(sum)

  # now do and print the HWE test
  hwe  <- HWE.test(data[,i])
  print(hwe)
}

