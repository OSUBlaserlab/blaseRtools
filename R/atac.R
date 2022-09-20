#' @title Split Out Peaks Data
#' @description Extracts Peaks data from 10X counts matrix and feature metadata and saves it as an alternate experiment in the assigned CDS.
#' @param cds CDS to split the Peaks data from
#' @return A cell data set
#' @rdname bb_split_atac
#' @export
#' @importFrom SingleCellExperiment splitAltExps swapAltExp
bb_split_atac <- function (cds) {
  check <- "Peaks" %in% bb_rowmeta(cds)$data_type
  stopifnot(`Your cds does not have any ATAC data.` = check)
  cds <-
    SingleCellExperiment::splitAltExps(cds,
                                       rowData(cds)$data_type,
                                       ref = "Peaks")
  cds <-
    SingleCellExperiment::swapAltExp(cds, name = "Gene Expression")
  cds <- as(object = cds, Class = "cell_data_set")
  return(cds)
}
