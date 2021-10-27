#' Convert windows filepath to linux-compatible.
#'
#' @description This will only work if the backslashes have already been escaped to double backslashes.  For example, readr automatically does this when reading in windows file paths from csv files. If you are using as a standalone function, run this in the terminal:  scan(what = "character, n = 1), paste the unmodified filepath at the blank line, and copy the modified file path as the argument, x.
#' @param x A character string filepath copied from Windows with escaped backslashes.
#' @return A linux-compatible filepath
#' @export
#' @import stringr
bb_fix_file_path <- function(x) {
  x <- str_replace_all(x,"\\\\","/")
  x <- str_replace(x,"X:/","~/network/X/")
  return(x)
}
