# $Id: plot.LD.R,v 1.3 2003/05/27 18:32:06 warnesgr Exp $

plot.LD.data.frame <- function(x,
                               digits=3,

                               sig=c(0,0.01, 0.025, 0.5, 0.1, 1),
                               celcol=heat.colors(length(sig)),
                               textcol="black",

                               marker,
                               which="D'",
                               distance,
                               ...)
  {
    oldpar <- par("mfrow")
    
    par(mfrow=c(1,2))

    LDtable(x, digits=digits, sig=sig, celcol=celcol, textcol=textcol, ...)
    LDplot(x, marker=marker, which=which, distance=distance, ...)
    
    par(mfrow=oldpar)
    invisible()
  }


LDtable <- function(x, 
                     sig=c(0,0.01, 0.025, 0.5, 0.1, 1),
                     celcol=heat.colors(length(sig)),
                     textcol="black",
                     digits=3,
                     ...)
  {
    dmat <- format( sign(x$"D") * x$"D'", digits=digits )
    tmp <- cut(x$"P-value", sig, include.lowest=TRUE)
    sigmat <- matrix(as.numeric(tmp),
                     nrow=nrow(x$"P-value"),
                     ncol=ncol(x$"P-value"))

    # remove blank row/column
    sigmat <- sigmat[-nrow(sigmat),-1]
    dmat <- dmat[-nrow(dmat),-1]
    n <- paste("(",x$n[-nrow(x$n),-1],")",sep="")
    p.value <- format.pval(x$"P-value"[-nrow(x$"P-value"),-1], digits=digits)
    
    nlev <- nlevels(tmp)

    image(x=1:ncol(sigmat), y=1:ncol(sigmat), z=sigmat[ncol(sigmat):1,],
          col=celcol, xlab="Marker 2\n\n", ylab="Marker 1",
          xaxt="n", yaxt="n",...)
    
    abline(v=-0.5 + 1:(ncol(sigmat)+1))
    abline(h=-0.5 + 1:(nrow(sigmat)+1))
    
    axis(3, 1:ncol(sigmat), colnames(dmat) )
    axis(2, 1:nrow(sigmat), rownames(dmat) )
    
    
    row <- matrix( 1:nrow(dmat), nrow=nrow(dmat), ncol=ncol(dmat), byrow=FALSE)
    col <- matrix( 1:nrow(dmat), nrow=nrow(dmat), ncol=ncol(dmat), byrow=TRUE)
        
    txtdat <- paste(dmat, n, p.value, sep="\n" )
    txtdat <- gsub("NA.*","",txtdat)
        
    text(x=nrow(dmat)-(row-1) ,
         y=col,
         txtdat,
         col=textcol,
         adj=c(0.5, 0.5)
         )
    
    text(x=1, y=1, "(-)D'\n(N)\nP-value", adj=c(0.5,0.5) )
    title(main="Linkage Disequilibrium\n")

    invisible()
  }



LDplot <- function(x, 
                   digits=3,
                   marker,
                   distance,
                   which="D'",
                   ...)
{
  if(missing(marker))
    marker <- colnames(x[[which]])
  else if (is.numeric(marker))
    marker <- colnames(x[[which]])[marker]
  
  datamat <- ifelse( is.na(x[[which]]), t(x[[which]]), x[[which]])

  if(which %in% c("D'","r") )
    diag(datamat) <- 1.0
  else if (which=="P-value")
    diag(datamat) <- 0.0
  
  dimnames(datamat) <- dimnames(x[[which]])
  
  if(missing(distance)) distance <- 1:ncol(datamat)
  distance <- matrix(distance, ncol=ncol(datamat), nrow=nrow(datamat),
                     byrow=TRUE)
  dimnames(distance) <- dimnames(datamat)
  
  matplot(x=t(distance[marker,,drop=FALSE]),
          t(datamat[marker,,drop=FALSE]),
          type="b", 
          xlab="Marker",
          ylab=paste("Linkage Disequilibrium: ", which, sep=""),
          xaxt="n",
          ... )
  
  axis(1, distance[1,], paste(1:ncol(datamat), colnames(datamat), sep=": " ))
  
  title("Pairwise Disequilibrium Plot")

  invisible()
}
