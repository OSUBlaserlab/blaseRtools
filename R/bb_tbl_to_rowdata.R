#' A Function To Add Tibble Columns To CDS RowData
#'
#' @param cds A cell data set object
#' @param min_tbl A tibble containing only the columns you want to add plus one column for joining.
#' @param join_col The column in min_tbl containing the join information for the cds rowData. Defaults to "feature_id".
#' @return A celldataset
#' @export
#' @import tidyverse SingleCellExperiment
bb_tbl_to_rowdata <- function (cds, min_tbl, join_col = "feature_id")
{
  cds_row_ids <- bb_rowmeta(cds = cds) %>%
    select(feature_id)
  min_tbl %>%
    rename(feature_id = contains(join_col))
  full_tbl <-
    left_join(cds_row_ids,
              min_tbl)
  stopifnot(
    "Joining min_tbl to rowData failed.  Check for duplicates in min_tbl." =
	    identical(full_tbl$feature_id, rownames(rowData(cds)))
  )
  cols_to_add <- full_tbl %>%
    select(-feature_id)
  for (i in 1:ncol(cols_to_add)) {
    rowData(cds)$new <- cols_to_add %>% pull(i)
    names(rowData(cds))[names(rowData(cds)) == "new"] <- colnames(full_tbl)[i]
  }
  return(cds)
}
