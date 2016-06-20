#' Examines if a package is already installed
#'
#' Examines if the package is already present
#' @export
#' @examples 
#' \\dontrun{
#' is_installed("GlobalVAR")}
#' }

is_installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
