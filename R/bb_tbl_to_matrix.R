#' Convert a wide-form tibble a matrix 
#' 
#' @param data A wide form tibble to convert to a matrix.  The first column will become the rownames. 
#' @return A matrix 
#' @export
#' @import dplyr
bb_tbl_to_matrix <- function(data) {
  data <- data %>%
    as.data.frame()
  rownames(data) <- data[,1]
  return(as.matrix(data[,-1]))
}
