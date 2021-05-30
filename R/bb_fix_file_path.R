#' Convert windows filepath to linux-compatible 
#' 
#' @param x A character string filepath copied from Windows
#' @return A linux-compatible filepath 
#' @export
#' @import stringr
#' @examples
bb_fix_file_path <- function(x) {
  x <- str_replace_all(x,"\\\\","/")
  x <- str_replace(x,"X:/","~/network/X/")
}
