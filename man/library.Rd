\name{library.pos}
\title{Loading and Listing of Packages (backported from 1.8.0)}
\alias{library.pos}
\description{
  \code{library.pos} loads add-on packages, allowing specification of
  position in the search path.
}
\usage{
library.pos(package, help, pos = 2, lib.loc = NULL, character.only = FALSE,
            logical.return = FALSE, warn.conflicts = TRUE,
            keep.source = getOption("keep.source.pkgs"),
            verbose = getOption("verbose"), version)
}
\arguments{
  \item{package, help}{\link{name} or character string giving the name
    of a package.}
  \item{pos}{the position on the search list at which to attach the
    loaded package. Note that \code{.First.lib} may attach other
    packages, and \code{pos} is computed \emph{after} \code{.First.lib}
    has been run.}
  \item{lib.loc}{a character vector describing the location of \R
    library trees to search through, or \code{NULL}.  The default value
    of \code{NULL} corresponds to all libraries currently known.}
  \item{character.only}{a logical indicating whether \code{package} or
    \code{help} can be assumed to be character strings.}
  \item{version}{A character string denoting a version number of the 
     package to be loaded.  If no version is given, a suitable default
     is chosen.}
  \item{logical.return}{logical.  If it is \code{TRUE},  \code{FALSE} or
      \code{TRUE} is returned to indicate success.}
  \item{warn.conflicts}{logical.  If \code{TRUE}, warnings are
    printed about \code{\link{conflicts}} from attaching the new
    package, unless that package contains an object \code{.conflicts.OK}.}
  \item{keep.source}{logical.  If \code{TRUE}, functions ``keep their
    source'' including comments, see argument \code{keep.source} to
    \code{\link{options}}.}
  \item{verbose}{a logical.  If \code{TRUE}, additional diagnostics are
    printed.}
}
\details{

  This function is backported from the R 1.8.0 development tree and is 
  only defined by gregmisc if this package is loaded in an earlier
  version of R.  In 1.8.0, this function is merely an alias for the
  standard \code{library} function.  The only important difference is
  the addition of the "pos" argument which allows an imported package to
  be loaded further down the search path than the top entry.

  See the documentation for \code{library} for more details.
}
\seealso{
  \code{\link[base]{library}}
}
\keyword{data}