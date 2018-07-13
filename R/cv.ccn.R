
cv.ccn <-
function(x, y, tom,KK , lambda = NULL, pred.loss = c("misclass", 
    "loss"),nfolds = 5, foldid, delta = 2, ...) {
    if (missing(pred.loss)) 
        pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
      ccn.object<-ccn(x, y,tom=tom, KK = KK , lambda= lambda,delta=delta, ...)
                                
         
    lambda <- gcdclust.object$lambda
    # predict -> coef
    nz <- sapply(coef(ccn.object, type = "nonzero"), length)
    if (missing(foldid)) 
        foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist <- as.list(seq(nfolds))
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
        which <- foldid == i
        y_sub <- y[!which]
        outlist[[i]] <-  ccn(x = x[!which, , drop = FALSE], 
                                  y = y_sub, tom=tom,KK = KK , lambda = lambda,delta=delta, ...)
        
    }
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(ccn.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, 
        pred.loss, delta))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
        cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, ccn.fit = ccn.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.ccn"
    obj
}
