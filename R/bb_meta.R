#' A Function to Return CDS Cell Metadata in Tibble Form
#'
#' @param cds A cell data set object
#' @param row_name Optional name to provide for cell unique identifier.  Defaults to "cell_id"
#' @return A tibble
#' @export
#' @import tidyverse SingleCellExperiment
bb_cellmeta <- function(cds, row_name = "cell_id") {
  SingleCellExperiment::colData(cds) %>%
    as_tibble(rownames = row_name)
}

#' A Function to Return CDS Row Metadata in Tibble Form
#'
#' @param cds A cell data set object
#' @param row_name Optional name to provide for cell unique identifier.  Defaults to "feature_id"
#' @return A tibble
#' @export
#' @import tidyverse SingleCellExperiment
bb_rowmeta <- function(cds, row_name = "feature_id") {
  SingleCellExperiment::rowData(cds) %>%
    as_tibble(rownames = row_name)
}
