nyblom <-function (y,x) {
# PURPOSE: computes the Nyblom (1989) test for parameter constancy against
# nonstationary alternatives (here lm) and also its
# heteroskedasticity-robust version (here rlm)
# From Gauss code of L. Vanessa Smith. See Dees, di Mauro, Pesaran, Smith (2007).  
# Converted to R by Alexandros Gabrielsen

yy = y
zz = x
k = NCOL(x)

# check if zz'*zz is posdef: if so avoid doing pseudoinverse
hh = tryCatch(chol(t(zz)%*%zz), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
if (length(hh)==1) {
  mzinv= ginv(t(zz)%*%zz)
} else {  
  xx = t(zz)%*%zz
  mzinv = diag(NROW(xx))%*%solve(xx)
  rm(xx)
}
rm(hh)
                                      
e=yy-zz%*%(mzinv%*%(t(zz)%*%yy))
seesq=(t(e)%*%e)/(NROW(e)-k)

ex=NULL
for (j in 1:NCOL(zz)) {
  ex_t = zz[,j]*e
  ex = cbind(ex, ex_t)
}
exs=apply(ex,2,cumsum)
       
v=as.numeric(seesq)*(t(zz)%*%zz)
rv=t(ex)%*%ex

hh = tryCatch(chol(v), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
if (length(hh)==1) {
 v1= ginv(v)
} else {  
  v1 = diag(NROW(v))%*%solve(v)
}
rm(hh)

hh = tryCatch(chol(rv), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
if (length(hh)==1) {
  v2= ginv(rv)
} else {  
  v2 = diag(NROW(rv))%*%solve(rv)
}
rm(hh)

lm=(v1%*%(t(exs)%*%exs))/NROW(e)
dglm=diag(lm)
lm=sum(dglm)

rlm=(v2%*%(t(exs)%*%exs))/NROW(e)
dgrlm=diag(rlm)
rlm=sum(dgrlm)

return = list(lm, rlm)
}
       
