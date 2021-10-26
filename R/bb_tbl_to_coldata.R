#' A Function To Add Cell-Level Metadata To CDS ColData
#'
#' @param cds A cell data set object
#' @param min_tbl A tibble containing only the columns you want to add plus one column for joining.  Must be cell-level metadata.  For sample level metadata add during loading with bb_load....
#' @param join_col The column in min_tbl containing the join information for the cds colData. Defaults to "cell_id".
#' @return A celldataset
#' @export
#' @import tidyverse SingleCellExperiment
bb_tbl_to_coldata <-
  function(cds, min_tbl, join_col = "cell_id") {
    full_tbl <-
      left_join(
        bb_cellmeta(cds),
        min_tbl,
        by = c("cell_id" = join_col)
      )
    stopifnot("Joining min_tbl to colData failed.  Check for duplicates in min_tbl." = identical(full_tbl$cell_id, rownames(colData(cds))))
    full_tbl <- full_tbl %>% select(-cell_id)
    for (i in 1:ncol(full_tbl)) {
      colData(cds)$new <- full_tbl %>% pull(i)
      names(colData(cds))[names(colData(cds)) == "new"] <-
        colnames(full_tbl)[i]
    }
    return(cds)
  }
