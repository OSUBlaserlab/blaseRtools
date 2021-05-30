#' A function to add tbl columns to cds colData
#'
#' @param cds A cell data set object
#' @param min_tbl A tibble containing only the columns you want to add plus one column for joining.
#' @param join_col The column containing the join information for the cds colData. Defaults to "barcode_sample".
#' @return A celldataset
#' @export
#' @import tidyverse SingleCellExperiment
#' @examples
bb_tbl_to_coldata <-
  function(cds, min_tbl, join_col = "barcode_sample") {
    full_tbl <-
      left_join(
        colData(cds) %>% as_tibble(rownames = "barcode_sample") %>% select("barcode_sample"),
        min_tbl,
        by = c("barcode_sample" = join_col)
      ) %>%
      select(-barcode_sample)
    for (i in 1:ncol(full_tbl)) {
      colData(cds)$new <- full_tbl %>% pull(i)
      names(colData(cds))[names(colData(cds)) == "new"] <-
        colnames(full_tbl)[i]
    }
    return(cds)
  }
