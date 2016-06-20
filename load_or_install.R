#' Check if a packaged in loaded or instaleed
#'
#' Finds whether to Load or Install a missing package. In case the package is missing it will install it from the CRAN.
#' @export
#' @examples 
#' \\dontrun{
#' }
load_or_install<-function(package_names){  
  for(package_name in package_names)  {    
    if(!is_installed(package_name)) {install.packages(package_name,repos="http://lib.stat.cmu.edu/R/CRAN")}
    library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE)  
  }
}
