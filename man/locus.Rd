% $Id: locus.Rd,v 1.12 2004/09/21 21:45:45 warnes Exp $

\name{locus}
\alias{locus}
\alias{gene}
\alias{marker}
\alias{is.gene}
\alias{is.locus}
\alias{is.marker}
\alias{print.gene}
\alias{print.locus}
\alias{print.marker}
\alias{as.character.locus}
\alias{as.character.gene}
\alias{as.character.marker}
\alias{getlocus}
\alias{getmarker}
\alias{getgene}
\alias{locus<-}
\alias{marker<-}
\alias{gene<-}


%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Create and Manipulate Locus, Gene, and Marker Objects}
\description{
  \code{locus}, \code{gene}, and \code{marker} create objects to store
  information, respectively, about genetic loci, genes, and markers.

  \code{is.locus}, \code{is.gene}, and \code{ismarker} test whether an
  object is a member of the respective class.

  \code{as.character.locus}, \code{as.character.gene},
  \code{as.character.marker} return a character string containing a
  compact encoding the object.
  
  \code{getlocus}, \code{getgene}, \code{getmarker} extract locus data
  (if present) from another object.

  \code{locus<-}, \code{marker<-}, and \code{gene<-} adds locus data to
  an object.

}
\usage{
  locus(name, chromosome, arm=c("p", "q", "long", "short", NA),
        index.start, index.end=NULL)

  gene(name, chromosome, arm=c("p", "q", "long", "short"),
       index.start, index.end=NULL)

  marker(name, type, locus.name, bp.start, bp.end = NULL,
         relative.to = NULL, ...)

  is.locus(x)

  is.gene(x)

  is.marker(x)

  as.character.locus(x, ...)

  as.character.gene(x, ...)

  as.character.marker(x, ...)

  getlocus(x, ...)

  locus(x) <- value

  marker(x) <- value

  gene(x) <- value

}
\arguments{
  \item{name}{character string giving locus, gene, or marker name}
  \item{chromosome}{integer specifying chromosome number (1:23 for humans).}
  \item{arm}{character indicating long or short arm of the chromosome.
    Long is be specified by "long" or "p".  Short is specified by
    "short" or "q".}
  \item{index.start}{integer specifying location of start of locus or
    gene on the chromosome.  }
  \item{index.end}{optional integer specifying location of end of locus or
    gene on the chromosome.  }
  \item{type}{character string indicating marker type, e.g. "SNP"}
  \item{locus.name}{either a character string giving the name of the
    locus or gene (other details may be specified using \code{...}) or a
    \code{locus} or \code{gene} object.}
  \item{bp.start}{start location of marker, in base pairs}
  \item{bp.end}{end location of marker, in base pairs (optional)}
  \item{relative.to}{location (optional) from which \code{bp.start} and
    \code{bp.end} are calculated. }
  \item{...}{parameters for \code{locus} used to fill in additional
    details on the locus or gene within which the marker is located. }
  \item{x}{an object of class \code{locus}, \code{gene}, or
    \code{marker}, or (for \code{getlocus}, \code{locus<-},
    \code{marker<-}, and \code{gene<-}) an object that may contain a locus
    attribute or field, notably a \code{genotype} object.}
  \item{value}{\code{locus}, \code{marker}, or \code{gene} object}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  Object of class \code{locus} and \code{gene}are lists with the
  elements:
  \item{name}{character string giving locus, gene, or marker name}
  \item{chromosome}{integer specifying chromosome number (1:23 for humans).}
  \item{arm}{character indicating long or short arm of the chromosome.
    Long is be specified by "long" or "p".  Short is specified by
    "short" or "q".}
  \item{index.start}{integer specifying location of start of locus or
    gene on the chromosome.  }
  \item{index.end}{optional integer specifying location of end of locus or
    gene on the chromosome.  }

  Objects of class \code{marker} add the additional fields:
  \item{marker.name}{character string giving the name of the marker}
  \item{bp.start}{start location of marker, in base pairs}
  \item{bp.end}{end location of marker, in base pairs (optional)}
  \item{relative.to}{location (optional) from which \code{bp.start} and
    \code{bp.end} are calculated. }
    
}
%\references{ ~put references to the literature/web site here ~ }
\author{Gregory R. Warnes \email{Gregory\_R\_Warnes@groton.pfizer.com} }

\seealso{
  \code{\link{genotype}},
}

\examples{
ar2  <- gene("AR2",chromosome=7,arm="q",index.start=35)
ar2

par  <- locus(name="AR2 Psedogene", 
              chromosome=1, 
              arm="q",
              index.start=32,
              index.end=42)
par

c109t  <- marker(name="C-109T",
                 type="SNP",
                 locus.name="AR2",
                 chromosome=7, 
                 arm="q", 
                 index.start=35,
                 bp.start=-109,
                 relative.to="start of coding region")
c109t

c109t  <- marker(name="C-109T",
                 type="SNP",
                 locus=ar2,
                 bp.start=-109,
                 relative.to="start of coding region")
c109t




example.data   <- c("D/D","D/I","D/D","I/I","D/D",
                    "D/D","D/D","D/D","I/I","")
g1  <- genotype(example.data, locus=ar2)
g1

getlocus(g1)

summary(g1)
HWE.test(g1)

g2  <- genotype(example.data, locus=c109t)
summary(g2)

getlocus(g2)

heterozygote(g2)
homozygote(g1)

allele(g1,1)

carrier(g1,"I")

heterozygote(g2)
}
\keyword{ misc }
%%\keyword{genetics}%-- one or more ...
