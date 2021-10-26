#' A Function To Add Tibble Columns To CDS RowData
#'
#' @param cds A cell data set object
#' @param min_tbl A tibble containing only the columns you want to add plus one column for joining.
#' @param join_col The column in min_tbl containing the join information for the cds rowData. Defaults to "feature_id".
#' @return A celldataset
#' @export
#' @import tidyverse SingleCellExperiment
bb_tbl_to_rowdata <-
  function(cds, min_tbl, join_col = "feature_id") {
    full_tbl <-
      left_join(
        bb_rowmeta(cds),
        min_tbl,
        by = c("feature_id" = join_col)
      )
    stopifnot("Joining min_tbl to rowData failed.  Check for duplicates in min_tbl." = identical(full_tbl$feature_id, rownames(rowData(cds))))
    full_tbl <- full_tbl %>% select(-feature_id)
    for (i in 1:ncol(full_tbl)) {
      rowData(cds)$new <- full_tbl %>% pull(i)
      names(rowData(cds))[names(rowData(cds)) == "new"] <-
        colnames(full_tbl)[i]
    }
    return(cds)
  }
