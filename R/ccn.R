ccn <-
function(x, y, tom=tom,  KK,  nlambda = 25, method = c("ccnL1", 
    "ccnscad", "ccnmcp"), 
                     lambda.factor = ifelse(nobs < nvars, 0.01, 
                                            1e-4), lambda = NULL, lambda2 = 0, pf = rep(1, nvars), pf2 = rep(1, nvars), exclude, 
                     dfmax = nvars + 1, pmax = min(dfmax * 1.2, nvars), standardize = TRUE, 
                     eps = 1e-5, maxit = 700, delta = 2) {
  #################################################################################

  method <- match.arg(method)
  this.call <- match.call()
  y <- drop(y)
  x <- as.matrix(x)
  KK<- as.integer(KK)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  if (is.null(vnames)) 
    vnames <- paste("V", seq(nvars), sep = "")
  if (length(y) != nobs) 
    stop("x and y have different number of observations")
  #################################################################################
  #parameter setup
  if (length(pf) != nvars) 
    stop("The size of L1 penalty factor must be same as the number of input variables")
  if (length(pf2) != nvars) 
    stop("The size of L2 penalty factor must be same as the number of input variables")
  if (lambda2 < 0) 
    stop("lambda2 must be non-negative")
  maxit <- as.integer(maxit)
  lam2 <- as.double(lambda2)
  pf <- as.double(pf)
  pf2 <- as.double(pf2)
  isd <- as.integer(standardize)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
 
 
  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) 
      stop("Some excluded variables out of range")
    jd <- as.integer(c(length(jd), jd))
  } else jd <- as.integer(0)

  #################################################################################
  #lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) 
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    #flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) 
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  #################################################################################
  if(missing(tom)){
    fit <- switch(method, 
                   ccnL1 = ccnclusterpathc(x, y, KK , nlam, flmin, 
                                            ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, 
                                            nobs, nvars, vnames),
                   ccnscad = ccnscadpathc(x, y, KK,nlam,flmin,ulam,isd,eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, 
                                                  nobs, nvars, vnames),
                   ccnmcp = ccnmcppathc(x, y, KK,nlam,flmin,ulam,isd,eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, 
                                     nobs, nvars, vnames)
                   
    )   
  }  
 else {  
   tom <- as.matrix(tom)
  fit <- switch(method, 
                 ccnL1 = hsvmclusterpathc(x, y, tom,KK , nlam, flmin, 
                                   ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, 
                                   nobs, nvars, vnames),
                 ccnscad = hsvmclusterscadpathc(x, y, tom, KK,nlam,flmin,ulam,isd,eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, 
                                   nobs, nvars, vnames),
                 ccnmcp = mcppathc(x, y, tom,KK,nlam,flmin,ulam,isd,eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, 
                                 nobs, nvars, vnames)
	
 ) 
 }                
  
  if (is.null(lambda)) 
    fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  #################################################################################
  class(fit) <- c("ccn", class(fit))
  fit
}
