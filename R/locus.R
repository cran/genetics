# $Id: locus.R,v 1.4 2002/11/12 05:31:20 warnesgr Exp $
#
# $Log: locus.R,v $
# Revision 1.4  2002/11/12 05:31:20  warnesgr
# - Fix mismatches between documentation and code that we generating
#   warning messages.
#
# Revision 1.3  2001/05/30 22:12:34  warnes
# Updated documentation mostly.  Added as.character.genotype().
#
# Revision 1.2  2001/05/22 18:58:10  warnes
# - Renamed "location.X" to "as.character.X" to better represent thier
#   function.  (Note as.character isn't generic at the moment, so these
#   functions must be called directly.)
#
# - added function getlocus() to grab the appropriate locus field or
#   attribute from an arbitrary object.
#
# Revision 1.1  2001/05/22 18:21:00  warnes
# Created classes for "locus", "gene", and "marker" along with basic
# methods is.X, print.X, and generic utility function "location", which
# formats the locus infromation for printing.
#
#

#                                        
#  > -----Original Message-----
#  > From: Friedrich Leisch [mailto:Friedrich.Leisch@ci.tuwien.ac.at]
#  > Sent: Wednesday, May 16, 2001 8:13 AM
#  > To: Gregory_R_Warnes@groton.pfizer.com; Georg Kainradl
#  > Subject: gene vs. marker
#  > 
#  > 
#  > 
#  > Hi Greg,
#  > 
#  > maybe I'm now making a complete fool of myself, 
	
#Not at all...  I just called my local geneticist to come up with the answers below. ;^)
#
#  > but what is 
#  > the exact
#  > difference between begin and end of your classes gene and marker? I
#  > mean for a marker m you define
#  > 
#  > m$begin
#  > m$end
#  > m$gene$begin
#  > m$gene$end
#  > 
#  > 
#  > and gene is a compulsory argument. I've only dealt with 
#  > marker data so
#  > far, so I simply don't know what your intention for the 
#  > gene slot is.
#
#Well, I was thinking of a hirearchical structure where markers lie within specific genes.  I'm learning that markers don't have to lie within known genes, though.  
#
#  > 
#  > Could we make "marker" inherit from "gene"?
#
#Hmm. Perhaps both should inherit from "locus" , which is more general than "gene";
#
#"locus" should have attributes: name, chromosome, arm, and index.start, 
#with optional attribute: index.end
#
#For example the gene AR2 is located at 7q35:
#   name="AR2"
#   chromosome=7, 
#   arm="q", 
#   index.start=35
#
#and a related pseudo-gene is located at 1q32-42:
#   name="AR2 Psedogene", 
#   chromosome=1, 
#   arm="q",
#   index.start=32,
#   index.end=42.
#
#"gene" would inherit from locus, but may add other information fields, say database references, etc.  (I'm not too clear on these.)
#
#Markers always relate to a specific locus, so perhaps it makes sense to inherit directly from locus.  "marker" would add required fields: type, bp.start (base-pair) 
#and have optional fields: bp.end, relative.to
#
#So, for a SNP we've been using here on AR2:
#   name="C-109T"
#   chromosome=7, 
#   arm="q", 
#   index.start=35
#   bp.start=-109
#   relative.to="start of coding region"
#
#One thing I'm a bit uncomfortable with is how to deal with names.  EG, this marker has a marker name "C-106T", but is located on a gene named "AR2".  How do we want to deal with the names.  Should we just have two fields, in "locus" the field name is "name", and in "marker" the name is "marker.name"? 
#
#  > 
#  > Reason why I'm asking is that we want to add an attribute "gene" or
#  > "marker" (using your class definitions) to class 
#  > "genotype", such that
#  > we can store the information of a linkage .par file in a 
#  > dataframe of
#  > genotypes. I think this is the most natural representation, 
#  > because we
#  > get subsetting etc. for free.
#  > 
#  > If "marker" inherits from "gene", we would need only one additional
#  > attribute in genotype, which could be either a marker or a gene.
#
#Makes sense, just call the slot "locus", and it can hold either.
#
#  > 
#  > I think we also need a chromosome slot in both marker and gene such
#  > that a list of markers (or genes) can be plotted directly.
#
#Yes.  Do you want to have separate a chromosome object, or just include the details directly in the structure as above?  The chromosome information is simple enough that it seems reasonable to just include it.
#
#-Greg
#
#  > 
#  > .f
#  > 
#

getlocus  <- function(x,...)
{
  if(is.locus(x))
    return(x)
  else if(!is.null(x$locus))
        return(x$locus)
  else if(!is.null(attr(x,"locus")))
       return(attr(x,"locus"))
  else
    NULL
}

locus  <- function(name, chromosome, arm=c("p","q","long","short"),
                   index.start, index.end=NULL)
  {
    if(missing(arm)) stop("Argument \"arm\" is missing, with no default")
    arm  <- match.arg( arm )
    
    object  <-  list()
    object$name  <- name
    object$chromosome <- chromosome
    
    
    object$arm  <- switch( arm, p="p", q="q", long="p", short="q")
    object$index.start  <- index.start
    object$index.end  <- index.end
    class(object)  <- "locus"
    return(object)
  }


gene  <-  function(name, chromosome, arm=c("p","q","long","short"),
                   index.start, index.end=NULL)
{
  object  <- locus(name, chromosome, arm, index.start, index.end)
  class(object)  <- c("gene","locus")
  object
}


marker <- function(name, type,
                   locus.name, bp.start, bp.end=NULL, relative.to=NULL,
                   ...
                   )
{
  if(is.locus(locus.name))
      object <- locus.name
  else
    object  <-  locus(locus.name, ...)
  
  object$marker.name  <- name
  object$type  <- type
  object$bp.start  <- bp.start
  object$bp.end  <- bp.end
  object$relative.to  <- relative.to
  class(object)  <- c("marker","locus")
  object
}

is.locus  <- function(x)
    inherits(x, "locus")

is.gene  <- function(x)
    inherits(x, "gene")

is.marker  <- function(x)
    inherits(x, "marker")



as.character.locus  <- function(x,...)
  {
    loc <- paste( x$chromosome, x$arm, x$index.start, sep="" )
    if( !is.null(x$index.end ) && x$index.start != x$index.end )
      loc  <- paste(loc, "-", x$index.end, sep="")
    loc
  }

as.character.gene  <- function(x,...)
  as.character.locus(x,...)

as.character.marker  <- function(x,...)
  {
    loc  <- as.character.locus(x)
    loc  <- paste(loc, ":", x$bp.start, sep="")
    if(!is.null(x$bp.end)) loc  <-  paste(loc, "-", x$bp.end, sep="")
    loc
  }

print.locus  <-  function(x,...)
  {
    cat("Locus: ", x$name, " (", as.character.locus(x), ")\n", sep="" )
  }

print.gene  <-  function(x,...)
  {
    cat("Gene: ", x$name, " (", as.character.locus(x), ")\n", sep="" )
  }

print.marker  <- function(x,...)
  {
    cat("Marker: ", paste(x$name,":",x$marker.name,sep=""),
        " (", as.character.marker(x), ")\tType: ",x$type,"\n", sep="" )
  }

