schow <- function(y,x,chow_cut) {
# PURPOSE: do Sequential Chow-tests for breaks. Three tests are constructed
#   (1) Quandt LR = SUP(F)
#   (2) Mean (F)
#   (3) Andrews-Ploberger = ln of mean of (exp(«F))
# and the corresponding heteroskedasticity-robust versions of all three
#--------------------------------------------------------------------------
#    Input:
#
#     y = lhv data
#     x = rhv data
#     ccut=endpoints for sequential chow regressions
#     btstrp=indicates whether bootstrap will be conducted
#--------------------------------------------------------------------------
#    Output:
#     suplr1  = sup(f)
#     meanlr1 = mean(f)
#     aplr1   = andrews/ploberger statistic
#     Hetero Robust Versions of these statistics
#     suplr2
#     meanlr2
#     aplr2
#    maxobs  = observation number corresponding to suplr
#--------------------------------------------------------------------------
#  % From Gauss code of L. Vanessa Smith. See Dees, di Mauro, Pesaran, Smith (2007).  
# Converted to R by Alexandros Gabrielsen
#**************************************************************************
  
y=as.matrix(y)
nobs=NROW(y)
ktrim=floor(chow_cut*nobs)
n1t=ktrim
n2t=nobs-ktrim
lr=array(0,c((n2t-n1t+1),2))

# full sample  Moment Matrices @
xy=t(x)%*%y
xx=t(x)%*%x
yy=t(y)%*%y


# check if zz'*zz is posdef: if so avoid doing pseudoinverse
hh = tryCatch(chol(xx), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
#if (!is.logical( hh )) {
if (length(hh)==1) {
  xxi= ginv(xx)
} else {  
  xxi = diag(NROW(xx))%*%solve(xx)
}
rm(hh)


e0=y-x%*%xxi%*%xy
ss0=t(e0)%*%e0
xe=NULL
for (j in 1:NCOL(x)) {
  xe_t = x[,j]*e0
  xe = cbind(xe, xe_t)
}
xxee=t(xe)%*%xe

x1y=t(x[1:(n1t-1),])%*%y[1:(n1t-1),]
x1x1=t(x[1:(n1t-1),])%*%x[1:(n1t-1),]
xxee1=t(xe[1:(n1t-1),])%*%xe[1:(n1t-1),]
yy1=t(y[1:(n1t-1),])%*%y[1:(n1t-1),]



i=n1t
while (i<= n2t) {

  x1y=x1y+t(t(x[i,])*y[i,])
  x1x1=x1x1+x[i,]%*%t(x[i,])
  xxee1=xxee1+xe[i,]%*%t(xe[i,])
  yy1=yy1+y[i,]%*%t(y[i,])
  x2x2=xx-x1x1
  x2y=xy-x1y
  yy2=yy-yy1
  xxee2=xxee-xxee1
  
  hh = tryCatch(chol(x1x1), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
  if (length(hh)==1) {
    x1x1i= ginv(x1x1)
  } else {  
   x1x1i = diag(NROW(x1x1))%*%solve(x1x1)
  }
  rm(hh)


  hh = tryCatch(chol(x2x2), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
  if (length(hh)==1) {
    x2x2i = ginv(x2x2)
  } else {  
    x2x2i = diag(NROW(x2x2))%*%solve(x2x2)
  }
  rm(hh)

  b1=t(x1x1i)%*%x1y
  b2=t(x2x2i)%*%x2y

  ssb=(yy1-t(x1y)%*%b1) + (yy2-t(x2y)%*%b2)
  lr[(i-n1t+1),1]=(nobs-2*NCOL(x))*((ss0-ssb)/ssb)
                     
# @ -- Hetero- Robust F-test -- @
                     
  v=x1x1i%*%(xxee1)%*%x1x1i + x2x2i%*%(xxee2)%*%x2x2i
                     
  hh = tryCatch(chol(v), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
  if (length(hh)==1) {
    v1 = ginv(v)
  } else {  
    v1 = diag(NROW(v))%*%solve(v)
  }
  rm(hh)

  lr[(i-n1t+1),2]=t(b1-b2)%*%v1%*%(b1-b2)
  i=i+1
}

lr1=lr[,1]
lr2=lr[,2]

# -- Non Robust -- 
suplr1=max(lr1)
maxobs=n1t+which(lr1==max(lr1))-1
meanlr1=mean(lr1)
aplr1=log(mean(exp(0.5*lr1)))

#  -- Hetero Robust -- 
suplr2=max(lr2)
meanlr2=mean(lr2)
aplr2=log(mean(exp(0.5*lr2)))

return = c(suplr1,meanlr1,aplr1,suplr2,meanlr2,aplr2,maxobs)
}
