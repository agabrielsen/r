neweywest <- function (X,res) {

# PURPOSE: computes Newey-West Heteroskedasticity and Autocorrelation
# Consistent standard errors
# From Gauss code of L. Vanessa Smith. See Dees, di Mauro, Pesaran, Smith (2007).  

    
s = NCOL(X)
T = NROW(res)
q = floor(4*(T/100)^(2/9))  # following Newey-West suggestion

for (v in (1:q)) {
  cw=array(0, c(s,s))
  iw=v+1
  while (iw<=NROW(X)) {
    qi=as.matrix(X[iw,])%*%res[iw,1]%*%res[(iw-v),1]%*%X[iw-v,]
    cw=cw+qi
    iw=iw+1
  }
  cw_new=(1-(v/(q+1)))*(cw+t(cw))

  if (v==1) {
    cw_new1=cw_new
  } else if (v>1) {
    cw_new1=cw_new1+cw_new
  }
}

lalagsbda= diag(diag(res%*%t(res)))

invXX = diag(NROW(t(X)%*%X))%*%solve(t(X)%*%X)

varcovNW= (T/(T-s))*(invXX%*%t(X)%*%lalagsbda%*%X%*%invXX +invXX%*%cw_new1%*%invXX) 

seNW=sqrt(diag(varcovNW))

return=list(seNW,q)
}
