gamm_cdf  <- function (x, dof) {
#' Gamma Distribution CDF
#'
#' Calculates the CDF of the Gamma Distribution
#' @param x an Nx1 vector
#' @param dof scalar describing the degress of freedom
#' @return gamm_cdf
#' @export
#' @author Alexandros Gabrielsen

if (dof<=0) {
  stop('degres-of-freedom misspecified')
}

cdf = pgamma(x,dof, lower.tail=TRUE, log.p=FALSE)
cdf[which(x<0)] = 0

return(cdf)
}
