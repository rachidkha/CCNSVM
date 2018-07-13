
standardize<- function(x,isd,ju)
{
x<-as.matrix(x)
dim<- dim(x)
nobs<- dim[1]
nvars<-dim[2]
isd<-as.integer(isd)
ju<-as.integer(ju)
 fit<- .Fortran("standard",nobs,nvars,as.double(x),mat=double(nobs * nvars),ju,isd,xmean=double(nvars),xnorm=double(nvars),maj=double(nvars)) 
 mat=matrix(fit$mat,nobs,nvars)
}
