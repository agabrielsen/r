kraplob <-function(y,x) {
#  PURPOSE: computes the Ploberger and Kramer's (1992) maximal OLS
# cumulative sum (CUSUM) statistic, here callede m1, and its mean-square 
# variant, here called m2.
# From Gauss code of L. Vanessa Smith. See Dees, di Mauro, Pesaran, Smith (2007).  

hh = tryCatch(chol(t(x)%*%x), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
#if (!is.logical( hh )) {
if (length(hh)==1) {
  invxx = ginv(t(x)%*%x)
  } else {  
    zz = t(x)%*%x
    invxx = diag(NROW(zz))%*%solve(zz)
    rm(zz)
  }

ehat=y - x%*%invxx%*%(t(x)%*%as.matrix(y))
vcv= t(ehat)%*%ehat/(NROW(y)-NCOL(x))
t1=cumsum(ehat/as.numeric(sqrt(NROW(ehat)*vcv)))
m1=max(abs(t1))
t2=t1*t1 
m2=mean(t2)
                                
return = c(m1, m2)
}
