#' Examines if a package is already installed
#'
#' @export
#' @examples 
#' \\dontrun{
#' }

is_installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
