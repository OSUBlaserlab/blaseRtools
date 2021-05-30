#' Remove rows that have duplicates in a given column 
#' 
#' @param data A tibble. 
#' @param column A column to deduplicate 
#' @return A deduplicated tibble
#' @export 
#' @import dplyr
#' @examples
bb_remove_dupes <- function(data,
			    column) {
  data <- data %>%
    group_by(!!sym(column)) %>%
    mutate(duplicate_flag = n() > 1) %>%
    filter(!duplicate_flag) %>%   
    select(-duplicate_flag)
  return(data)
}
