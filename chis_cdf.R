#' CDF of Chi-Squared Distribution
#'
#' CDF of Chi-Squared Distribution
#' @param x a Nx1 vector
#' @param dof degress of freedom
#' @return chis_cdf
#' @export
#' @author Alexandros Gabrielsen <agabrielsen01@@outlook.com>

chis_cdf <- function(x, dof) {

if (dof<=0) {
  stop('degres-of-freedom misspecified')
}
cdf = gamm_cdf(x/2,dof*0.5) # pgamma(x/2,dof*0.5)
return(cdf)
}
