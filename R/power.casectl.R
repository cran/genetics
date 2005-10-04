# power calculation for case-control design (default: case:control = 1:1)
# Author: Michael Man
#   Date: May 5, 2004
#      N: total number of subjects
#  gamma: relative risk in multiplicative model;
#         not used in Dominant or Recessive model (assume A as protective allele)
#      p: frequency of A allele
#     kp: prevalence of disease
#  alpha: significance level
#     fc: fraction of cases
#     pi: probability of 'aa' genotype has the disease
#   minh: mode of inheritance
# reference: Long, A. D. and C. H. Langley (1997). Genetic analysis of complex traits. Science 275: 1328.
#            Agresti, A. (2002) Categorical Data Analysis. Second Edition, p243.
# ( modified from pbsize{gap} )
# requirement: It is recommended to use R 1.9.0 or above.  The function 'qchisq' in earlier version
#              has problem with large noncentrality parameter.


# under HWE                       AA       Aa       aa       
#   fHW = p(genotype)        = c( p^2,     2pq,     q^2 )     

# model specification 
#   f.mod = relative risk    = c(gamma^2,  gamma,    1 )    # multiplicative model
#   f.mod =                  = c(  0,        0,      1 )    #       dominant model
#   f.mod =                  = c(  0,        1,      1 )    #      recessive model

# conditional prob.
#   p(D|genotype) = f.mod*pi = c(gamma^2,  gamma,    1 )*pi

# population joint prob. (f.mod = 1 under Ho)
#   Case     p(D,     genotype) = p(genotype)*     p(D|genotype)  = fHW*   f.mod*pi
#   Control  p(D_not, genotype) = p(genotype)*(1 - p(D|genotype)) = fHW*(1-f.mod*pi)

# population conditional prob. (f.mod = 1 under Ho)
#   Case     p(genotype|D)     = p(D    , genotype)/P(D    ) = P(D    , genotype)/sum(P(D    , genotype)) = fHW*   f.mod*pi  /    sum(fHW*f.mod*pi)
#   Control  p(genotype|D_not) = p(D_not, genotype)/P(D_not) = P(D_not, genotype)/sum(P(D_not, genotype)) = fHW*(1-f.mod*pi) / (1-sum(fHW*f.mod*pi))

# sample or allocation probability
#   1:1 case-control design  p(D|Sample) = fc = 1/2
#   1:2 case-control design                fc = 1/3
#   a prospective design                   fc = sum(fHW*f.mod*pi)

# sample joint prob. (f.mod = 1 under Ho)
# for prospective design, this is the same as population joint prob. since 'fc' cancels out with 'sum(fHW*f.mod*pi)' 
#   Case     p(genotype,D    |sample) = p(genotype|D    )*     p(D|Sample)  =    fc *fHW*   f.mod*pi  /    sum(fHW*f.mod*pi)
#   Control  p(genotype,D_not|sample) = p(genotype|D_not)*(1 - p(D|Sample)) = (1-fc)*fHW*(1-f.mod*pi) / (1-sum(fHW*f.mod*pi))


power.casectrl <- function (N, gamma = 4.5, p = 0.15, kp=0.1, alpha=.05, fc=0.5,
                             minh=c('multiplicative', 'dominant', 'recessive')) 
{
  minh <- match.arg(minh)
  if ( !all(gamma > 0, N > 0) ) stop('N and gamma must be greater than 0')
  if ( min(p, kp, alpha, fc) <= 0 | max(p, kp, alpha, fc) >=1 ) stop('p, kp, alpha, and fc must be between 0 and 1.') 
  f.mod <- switch(minh,
         multiplicative = c(gamma^2, gamma, 1),
         dominant       = c(      0,     0, 1),
         recessive      = c(      0,     1, 1)  ) 
  q <- 1 - p
  fhw <- c(p^2, 2*p*q, q^2)
  pi <- kp/sum(f.mod*fhw)
  if (pi <= 0 | pi >=1) {
    warning('The combination of p, kp, and gamma produces an unrealistic value of pi.')
    ret <- NA
  } else {
    fe  <- rbind(fhw, fhw)
    dimnames(fe) <- list(c("Case", "Control"), c("AA", "Aa", "aa"))
    f <- fe*rbind(f.mod*pi, 1-f.mod*pi)
    Pct <- apply(f, 1, sum)
    f2 <-  f *c(fc, 1-fc)/Pct   # normalize the frequencies for each row
    fe2 <- fe*c(fc, 1-fc)  
    fe2; apply(fe2, 1, sum); f2; apply(f2, 1, sum)  
    lambda <- sum((f2-fe2)^2/fe2)*N
    ret <- 1 - pchisq(qchisq(1-alpha, df=1), df=1, ncp=lambda, lower.tail=T)
  }
  ret
}


power.casectrl.plot <- function (N, gamma=1.6, p=1:9/10, kp=0.1, alpha=0.05, fc=0.5,
                                 minh=c('multiplicative', 'dominant','recessive'),
                                 Nsnp=1, vary=c('prevalence','SNPs'), ylim=c(0,1), PLOT=T, ... )
{
  minh <- match.arg(minh)
  vary <- match.arg(vary)
  if (length(p)<2) stop('Must have more than 1 value in p.')
  if (length(kp) > 1 & length(Nsnp) > 1) stop("Nsnps and kp can't be all > 1.")
  if (vary=='prevalence') {
    cmd <- expression(tapply(p, p, function(x, ...) power.casectrl(p=x,...), N=N, gamma=gamma, kp=kp[j], alpha=alpha/Nsnp, fc=fc, minh=minh))
    Xvary <- kp 
  } else if (vary=='SNPs') {
    cmd <- expression(tapply(p, p, function(x, ...) power.casectrl(p=x,...), N=N, gamma=gamma, kp=kp, alpha=alpha/Nsnp[j], fc=fc, minh=minh))
    Xvary <- Nsnp 
  }
  J <- length(Xvary)
  ret <- matrix(NA, nc=J, nr=length(p))
  colnames(ret) <- paste(vary, '=', Xvary)
  for (j in 1:J) ret[,j] <- eval(cmd)

  if (PLOT) {
    nc <- 1:ncol(ret)
    subt <- paste("( RR", gamma, "; total subjects", N,"; SNPs", Nsnp[1], "; prevalence", kp[1],
                   "; mode of inheritance:", minh, "; overall sig.level", alpha, ")" )
    matplot(p, ret, type="l", ylim=ylim, lty=1, col=nc, xlab="Allele Frequency", ylab="Power", sub=subt, ...)
    abline(h=c(.8), lty=1)
    legend( locator(1), colnames(ret), lty=1, col=nc )
  }
  ret
}



### power calculation for case-control design (case:control = 1:1)
### old version, less flexible
power.casectrl.old <- function (N, gamma = 4.5, p = 0.15, kp=.1, alpha=.05) # modified from pbsize{gap}
{
  q <- 1 - p
  pi <- kp/(gamma * p + q)^2
  f <- array(99, dim=c(2,3))
  dimnames(f) <- list(c("Case", "Control"), c("AA", "Aa", "aa"))
  f['Case', 'AA'] <-   pi*(gamma*p)^2
  f['Case', 'Aa'] <- 2*pi* gamma*p *q
  f['Case', 'aa'] <-   pi* q^2
  f['Control', 'AA'] <-   p^2
  f['Control', 'Aa'] <- 2*p*q
  f['Control', 'aa'] <-   q^2

  f <- fe <- f/2
  fe['Case',] <- fe['Control',]
  Pct <- pi*(gamma*p+q)^2
  f['Control',] <- (f['Control',] - f['Case',])/(1-Pct)
  f['Case',] <- f['Case',]/Pct
  fe
  f
  sum(fe); sum(f[1,]); sum(f[2,]) 
  lambda <- sum((f-fe)^2/fe)*N
  ret <- pchisq(qchisq(1-alpha, df=1), df=1, ncp=lambda, lower.tail=F)
  ret
}

### simple simulation for two group design
pw <- function(n1, n2=n1*(1-fc)/fc, fc=.5, pi=0, me1=50, me2=45, sd1=10, sd2=10, TEST=F){
  covm <- matrix(c(1,    pi,    pi,    1   ), nr=2)*
          matrix(c(sd1^2, sd1*sd2, sd1*sd2, sd2^2), nr=2)
  x1 <- rmvnorm(n=n1, mean=c(me1,me1), sigma=covm)
  x2 <- rmvnorm(n=n2, mean=c(me1,me2), sigma=covm)
  x <- data.frame(rbind(x1,x2), Trt=c(rep(0,n1), rep(1,n2)))
  colnames(x) <- c('X', 'Y', 'Trt')
  mod <- lm(Y~X+Trt, data=x)
  if (TEST) {
    print( summary(mod) )
    plot(Y~X+Trt, data=x)
  }
  ret <- anova(mod, test='F')['Trt','Pr(>F)']
  ret
}

# power calculation for studies using baseline measure 
#   - simulation: continuous response, baseline and genotype as covariate (ANCOVA)
#   - can specify various modes of inheritance ('additive', 'dominant','recessive')
#   - use compound symmetry for covariance matrix 
# Author: Michael Man
#   Date: June 22, 2004
#      N: total number of subjects
#      p: frequency of A allele
#  alpha: significance level                     (used only in 'power.genotype.conti')
#    Rep: number of iterations to generate power (used only in 'power.genotype.conti')
#     pi: correlation coefficient
#    me1: mean of control group
#  delta: treatment/genotype effect
#  sd1/2: standard deviation of the control and treatment groups
#   minh: mode of inheritance
# genotype.delta: the effect due to individual genotype effect or overall effect
# Factor: whether treat 'Trt' as a factor in 'x'
#   TEST: debug
# reference:
#   Frison and Pocock (1992) "Repeated measures in clinical trials: analysis using mean summary statistics
#     and its implications for design" Statistics in Medicine 11:1685-1704
#   Vickers (2001) "The use of percentage change from baseline as an outcome in a controlled trial is
#     statistically inefficient: a simulation study" BMC Med Res Methodol. 2001; 1 (1): 6
# requirement: need library 'mvtnorm'

simu.genotype.conti <- function (N, p=0.15, pi=0, me1=50, me2=me1, delta=-5, sd1=10, sd2=10, TEST=F,
                                 minh=c('additive', 'dominant','recessive'), genotype.delta=T, Factor=F) 
{
  minh <- match.arg(minh)
  if ( min(N, sd1, sd2)<0 ) stop('N, sd1, and sd2 must be greater than 0')
  if ( p<=0 | abs(pi)<0 | p>=1 | abs(pi)>1 ) stop('p and abs(pi) must be between 0 and 1.') 
  f.mod <- switch(minh,
         dominant       = c(      0,     1, 1),
         additive       = c(      0 ,  0.5, 1),
         recessive      = c(      0,     0, 1)  ) 
  q <- 1 - p
  fhw <- c(q^2, 2*p*q, p^2)  # major allele first
  nhw <- round(N*fhw)
  if (sum(nhw)!=N) nhw[3] <- N-sum(nhw[1:2])
  covm <- matrix(c(1,    pi,    pi,    1   ), nr=2)*
          matrix(c(sd1^2, sd1*sd2, sd1*sd2, sd2^2), nr=2)
  if (!genotype.delta) delta <- delta/sum(fhw*f.mod)  # convert to overall delta due to all genotypes
  for (i in 1:3) {
    if (nhw[i]!=0) {
      assign(paste('x',i,sep=''), rmvnorm(n=nhw[i], mean=c(me1,me2+f.mod[i]*delta), sigma=covm))
      assign(paste('t',i,sep=''), rep(f.mod[i],nhw[i]))
    } else {
      assign(paste('x',i,sep=''), NULL)
      assign(paste('t',i,sep=''), NULL)
    }
  }
  x <- data.frame(rbind(x1,x2,x3), Trt=c(t1,t2,t3))
  if (Factor) x$Trt <- as.factor(x$Trt)
  colnames(x) <- c('X', 'Y', 'Trt')
  mod <- lm(Y~X+Trt, data=x)   # ANCOVA
  if (TEST) {
    print( summary(mod) )
    plot(Y~X+Trt, data=x)
  }
  ret <- anova(mod, test='F')['Trt','Pr(>F)']
  ret
}

### power calculation
power.genotype.conti <- function(N, Rep=100, alpha=.05, ...){
  pval <- sapply(rep(N, Rep), FUN=simu.genotype.conti, ...)
  power <- length(pval[pval<=alpha])/length(pval)
  power
}


