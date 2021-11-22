#' A Function To Add Cell-Level Metadata To CDS ColData
#'
#' @param cds A cell data set object
#' @param min_tbl A tibble containing only the columns you want to add plus one column for joining.  Must be cell-level metadata.  For sample level metadata add during loading with bb_load....
#' @param join_col The column in min_tbl containing the join information for the cds colData. Defaults to "cell_id".
#' @return A celldataset
#' @export
#' @import tidyverse SingleCellExperiment
bb_tbl_to_coldata <- function (cds, min_tbl, join_col = "cell_id")
{
  cds_cell_ids <- bb_cellmeta(cds = cds) %>%
    select(cell_id)
  min_tbl %>%
    rename(cell_id = contains(join_col))
  full_tbl <-
    left_join(cds_cell_ids,
              min_tbl)
  stopifnot(
    "Joining min_tbl to colData failed.  Check for duplicates in min_tbl." =
	    identical(full_tbl$cell_id, rownames(colData(cds)))
  )
  cols_to_add <- full_tbl %>%
    select(-cell_id)
  for (i in 1:ncol(cols_to_add)) {
    colData(cds)$new <- cols_to_add %>% pull(i)
    names(colData(cds))[names(colData(cds)) == "new"] <- colnames(cols_to_add)[i]
  }
  return(cds)
}
