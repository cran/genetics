\name{print.LD}
\alias{print.LD}
\alias{print.LD.data.frame}
\alias{plot.LD.data.frame}
\alias{LDtable}
\alias{LDplot}
\title{Textual and graphical display of linkage disequilibrium (LD) objects}
\description{
  Textual and graphical display of linkage disequilibrium (LD) objects
}
\usage{
print.LD(x, digits = getOption("digits"), ...)
print.LD.data.frame(x, digits=getOption("digits"),
                    which=c("D", "D'", "Corr.", "Chi-Square",
                            "P-Value", "N"), rowsep, show.all=FALSE,
                    ...)

plot.LD.data.frame(x,digits=3, sig=c(0,0.01, 0.025, 0.5, 0.1, 1),
                   celcol=heat.colors(length(sig)), textcol="black",
                   marker, which="D'", distance,  ...)


LDtable(x, sig=c(0,0.01, 0.025, 0.5, 0.1, 1),
        celcol=heat.colors(length(sig)), textcol="black",
        digits=3, ...)

LDplot(x, digits=3, marker=colnames(x$D), distance, which="D'", ...)
}
\arguments{
  \item{x}{LD or LD.data.frame object}
  \item{digits}{Number of significant digits to display}
  \item{which}{Name(s) of LD information items to be displayed.}
  \item{rowsep}{Separator between rows of data, use \code{NULL} for no separator.}
  \item{sig}{P-value cutoffs points for colorizing LDtable}
  \item{celcol}{Colors fore each P-value cutoff given in \code{sig} for
  LDtable}
  \item{textcol}{Color for text labels for LDtable}
  \item{marker}{Marker used as 'comparator' on LDplot.  If
    omitted separate lines for each marker will be displayed.}
  \item{distance}{Marker location, used for locating of markers on
    LDplot.}
  \item{show.all}{If TRUE, show all rows/columns of matrix. Otherwise
    omit completely blank rows/columns.}
  \item{\dots}{Optional arguments (\code{plot.LD.data.frame} passes
    these to \code{LDtable} and \code{LDplot})}
}
%\details{
%}
\value{
  None.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Gregory R. Warnes \email{gregory\_r\_warnes@groton.pfizer.com} }
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

# NB: more markers makes prettier plots!
print(ldt)                               # textual display
LDtable(ldt)                            # graphical color-coded table
LDplot(ldt, distance=c(124, 834, 927))  # LD plot vs distance

plot(ldt) # do LDtable & LDplot on in a single graphics window
}
\keyword{misc}

