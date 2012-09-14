% $Id: print.LD.Rd 1353 2012-08-14 14:23:17Z warnes $

\name{print.LD}
\alias{print.LD}
\alias{print.LD.data.frame}
\alias{summary.LD.data.frame}
\alias{print.summary.LD.data.frame}
\alias{plot.LD.data.frame}
\alias{LDtable}
\alias{LDplot}
\title{Textual and graphical display of linkage disequilibrium (LD) objects}
\description{
  Textual and graphical display of linkage disequilibrium (LD) objects
}
\usage{

\method{print}{LD}(x, digits = getOption("digits"), \dots)
\method{print}{LD.data.frame}(x, \dots)

\method{summary.LD}{data.frame}(object, digits = getOption("digits"),
                      which = c("D", "D'", "r", "X^2", "P-value", "n", " "),
                      rowsep, show.all = FALSE, \dots)
\method{print}{summary.LD.data.frame}(x, digits = getOption("digits"), \dots)

\method{plot}{LD.data.frame}(x,digits=3, colorcut=c(0,0.01, 0.025, 0.5, 0.1, 1),
                   colors=heat.colors(length(colorcut)), textcol="black",
                   marker, which="D'", distance,  \dots)


LDtable(x, colorcut=c(0,0.01, 0.025, 0.5, 0.1, 1),
        colors=heat.colors(length(colorcut)), textcol="black",
        digits=3, show.all=FALSE, which=c("D", "D'", "r", "X^2",
        "P-value", "n"), colorize="P-value", cex, \dots)

LDplot(x, digits=3, marker, distance, which=c("D", "D'", "r", "X^2",
       "P-value", "n", " "), \dots )
}
\arguments{
  \item{x,object}{LD or LD.data.frame object}
  \item{digits}{Number of significant digits to display}
  \item{which}{Name(s) of LD information items to be displayed}
  \item{rowsep}{Separator between rows of data, use \code{NULL} for no
    separator.}
  \item{colorcut}{P-value cutoffs points for colorizing LDtable}
  \item{colors}{Colors for each P-value cutoff given in \code{colorcut} for
    LDtable}
  \item{textcol}{Color for text labels for LDtable}
  \item{marker}{Marker used as 'comparator' on LDplot.  If
    omitted separate lines for each marker will be displayed}
  \item{distance}{Marker location, used for locating of markers on
    LDplot.}
  \item{show.all}{If TRUE, show all rows/columns of matrix. Otherwise
    omit completely blank rows/columns.}
  \item{colorize}{LD parameter used for determining table cell colors}
  \item{cex}{Scaling factor for table text. If absent, text will be
    scaled to fit within the table cells.}
  \item{\dots}{Optional arguments (\code{plot.LD.data.frame} passes
    these to \code{LDtable} and \code{LDplot})}
}
%\details{
%}
\value{
  None.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Gregory R. Warnes \email{greg@warnes.net} }
\seealso{ \code{LD}, \code{genotype}, \code{HWE.test} }
\examples{


g1 <- genotype( c('T/A',    NA, 'T/T',    NA, 'T/A',    NA, 'T/T', 'T/A',
                  'T/T', 'T/T', 'T/A', 'A/A', 'T/T', 'T/A', 'T/A', 'T/T',
                     NA, 'T/A', 'T/A',   NA) )

g2 <- genotype( c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
                  'C/A', 'C/C', 'C/A', 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
                  'C/A', 'C/A', 'C/A', 'A/A') )


g3 <- genotype( c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
                  'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A', 'T/T',
                  'T/A', 'T/A', 'T/A', 'T/T') )
data <- makeGenotypes(data.frame(g1,g2,g3))

# Compute & display  LD for one marker pair
ld <- LD(g1,g2)
print(ld)

# Compute LD table for all 3 genotypes
ldt <- LD(data)

# display the results
print(ldt)                               # textual display
LDtable(ldt)                            # graphical color-coded table
LDplot(ldt, distance=c(124, 834, 927))  # LD plot vs distance

# more markers makes prettier plots!
data <- list()
nobs <- 1000
ngene <- 20
s <- seq(0,1,length=ngene)
a1 <- a2 <- matrix("", nrow=nobs, ncol=ngene)
for(i in 1:length(s) )
{

  rallele <- function(p) sample( c("A","T"), 1, p=c(p, 1-p))

  if(i==1)
    {
      a1[,i] <- sample( c("A","T"), 1000, p=c(0.5,0.5), replace=TRUE)
      a2[,i] <- sample( c("A","T"), 1000, p=c(0.5,0.5), replace=TRUE)
    }
  else
    {
      p1 <- pmax( pmin( 0.25 + s[i] * as.numeric(a1[,i-1]=="A"),1 ), 0 )
      p2 <- pmax( pmin( 0.25 + s[i] * as.numeric(a2[,i-1]=="A"),1 ), 0 )
      a1[,i] <- sapply(p1, rallele )
      a2[,i] <- sapply(p2, rallele )
    }

  data[[paste("G",i,sep="")]] <- genotype(a1[,i],a2[,i])
}
data <- data.frame(data)
data <- makeGenotypes(data)

ldt <- LD(data)
plot(ldt, digits=2, marker=19) # do LDtable & LDplot on in a single
                               # graphics window
}
\keyword{misc}
