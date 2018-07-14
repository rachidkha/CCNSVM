
checkv<-function(x)
{
  x<-as.matrix(x)
  dim<- dim(x)
  nobs<- dim[1]
  nvars<-dim[2]
  fit<-.Fortran("chkvars",nobs,nvars,as.double(x),ju=integer(nvars))
  fit$ju
}
