# $Id: print.LD.R,v 1.4 2003/05/27 21:18:06 warnesgr Exp $

print.LD <- function(x, digits=getOption("digits"), ...)
  {
    saveopt <- options("digits")
    options(digits=digits)
    cat("\n")
    cat("Pairwise LD\n")
    cat("-----------\n")

    est <- t(as.matrix( c(D=x$"D","D'"=x$"D'","Corr"=x$"r")))
    rownames(est) <- "Estimates:"
    print(est)
    cat("\n")

    test <- t(as.matrix( c("X^2"=x$"X^2", "P-value"=x$"P-value",
                           "N"=x$"n") ) )
    rownames(test) <- "LD Test:"
    print(test)
    cat("\n")

    options(saveopt)
    invisible(x)
  }


print.LD.data.frame <- function(x, digits=getOption("digits"),
                                which=c("D", "D'", "r", "X^2",
                                        "P-value", "n", " "),
                                rowsep, show.all=FALSE,
                                ...)
  {

    if(missing(rowsep))
      if(length(which)==1)
        rowsep <- NULL
      else
        rowsep <- " "

    if(is.null(rowsep))
      blank <- NULL
    else
      blank <- matrix(rowsep, ncol=ncol(x$"D"), nrow=nrow(x$"D"))
    


    saveopt <- options("digits")
    options(digits=digits)

    
    pdat <- list()
    for(name in which)
      pdat[[name]] <- x[[name]]
    
    tab <- interleave(
                      D = pdat$"D",
                      "D'" = pdat$"D'",
                      "Corr." = pdat$"r",
                      "X^2"= pdat$"X^2",
                      "P-value" = pdat$"P-value",
                      "n" = pdat$"n",
                      " "=blank,
                      sep=" "
                      )

    statlist <- which[ ! (which %in% c("P-value", "n", " ") ) ]
    statlist[statlist=="X^2"] <- "X\\^2"

    formatlist <- sapply( statlist, function(x) grep(x, rownames(tab) ) )
    formatlist <- unique(sort(unlist(formatlist)))
    
    pvallist   <- grep( "P-value", rownames(tab) )
    
    tab[formatlist,] <- formatC(as.numeric(tab[formatlist,]), digits=digits,
                                format="f")
    tab[pvallist,] <- apply(x$"P-value", c(1,2),
                            function(x)trim(format.pval(x, digits=digits)))
    
    tab[trim(tab)=="NA"] <- NA

    if(!show.all)
      {
         # drop blank row/column
        entrylen <- nrow(tab)/nrow(x$n)
        tab <- tab[1:(nrow(tab) - entrylen),-1]
      }

    cat("\n")
    cat("Pairwise LD\n")
    cat("-----------\n")
    
    print.matrix(tab, digits=digits, quote=FALSE, na.print="    ", right=TRUE) 

    cat("\n")
    
    options(saveopt)
    invisible(tab)
  }
