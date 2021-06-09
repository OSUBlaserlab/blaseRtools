#' Rejoin qc and doubletfinder data to a cds object 
#' 
#' @param cds A cell data set object to rejoin to
#' @param qc_data A table of cell barcodes with qc data.  Can be extracted from bb_qc with purrr::map(qc_result, 1) 
#' @param doubletfinder_data The doubletfinder result tbl
#' @return A cell data set object with qc and doubletfinder data 
#' @export
#' @import tidyverse monocle3
bb_rejoin <- function(cds, qc_data, doubletfinder_data) {
  cds_tbl <- as_tibble(colData(cds))
  cds_tbl <- left_join(cds_tbl, qc_data)
  cds_tbl <- left_join(cds_tbl, doubletfinder_data)
  cds_df <- as.data.frame(cds_tbl)
  row.names(cds_df) <- cds_df$barcode
  cds <- new_cell_data_set(
    expression_data = cds@assays@data$counts,
    cell_metadata = cds_df,
    gene_metadata = rowData(cds)
  )
  
  
  return(cds)
}
