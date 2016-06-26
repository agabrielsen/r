intercept.model <- function(Y) {
# Fit an intercept model
# Author Alexandros Gabrielsen
# Fitting an intercept model
bb = tryCatch(arima(Y, c(0,0,0),method="ML", optim.method="BFGS",optim.control = list(maxit = 2500, reltol=1e-8)), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
if( !is.logical( bb ) ) {
  return = bb
} else { return =  NULL }
}
