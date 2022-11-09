#' @title Get Cell Metadata
#' @description Take a cell_data_set object or a Seurat object and return the cell metadata in the form of a tibble.  The unique cell identifier column is labeled cell_id by default.  Prior versions of this function would only accept a cell_data_set.  The input argument has been changed from cds to obj to reflect the fact that Seurat objects are now also accepted.
#' @param obj A cell_data_set or Seurat object.
#' @param row_name Optional name to provide for cell unique identifier, Default: 'cell_id'
#' @param cds Provided for compatibility with prior versions, Default: NULL
#' @return A tibble
#' @details If a value is supplied for cds, a warning will be issued and the function will pass the value of cds to obj.
#' @rdname bb_cellmeta
#' @export
#' @importFrom tibble as_tibble
bb_cellmeta <- function(obj, row_name = "cell_id", cds = NULL) {
  cds_warn(cds)
  obj_stop(obj)
  if ("cell_data_set" %in% class(obj))
    res <- tibble::as_tibble(x = obj@colData, rownames = row_name)
  if ("Seurat" %in% class(obj))
    res <- tibble::as_tibble(x = obj@meta.data, rownames = row_name)
  return(res)
}



#' @title Get Feature/Gene Metadata
#' @description Take a cell_data_set or Seurat object and return the gene/feature metadata in the form of a tibble.  RNA is used as the default assay.
#' @param obj A cell_data_set or Seurat object
#' @param row_name Optional name to provide for feature unique identifier, Default: 'feature_id'
#' @param experiment_type The experiment type to display.  Applies only to cds objects.  Commonly will be either "Gene Expression" or "Antibody Capture", Default: 'Gene Expression'
#' @param assay For a Seurat object, th feature assay to return.  CDS objects with alternative experiments are not supported, Default: 'RNA'
#' @param cds Provided for compatibility with prior versions, Default: NULL
#' @return At tibble.
#' @details If a value is supplied for cds, a warning will be issued and the function will pass the value of cds to obj.
#' @rdname bb_rowmeta
#' @export
#' @importFrom SingleCellExperiment mainExpName altExpNames swapAltExp rowData
#' @importFrom cli cli_abort
#' @importFrom tibble as_tibble
bb_rowmeta <- function(obj,
                       row_name = "feature_id",
                       experiment_type = "Gene Expression",
                       assay = "RNA",
                       cds = NULL) {
  cds_warn(cds)
  obj_stop(obj)
  if ("cell_data_set" %in% class(obj)) {
    # check to be sure experiment_type is available
    all_exps <- c(
      SingleCellExperiment::mainExpName(obj),
      SingleCellExperiment::altExpNames(obj)
    )
    if (experiment_type %notin% all_exps)
      cli::cli_abort("The requested experiment name is not available.")
    if (experiment_type != "Gene Expression")
      obj <-
        as(SingleCellExperiment::swapAltExp(obj, name = experiment_type),
           Class = "cell_data_set")


    res <-
      tibble::as_tibble(x = SingleCellExperiment::rowData(obj), rownames = row_name)
  }
  if ("Seurat" %in% class(obj))
    res <-
    tibble::as_tibble(x = obj[[assay]][[]], rownames = row_name)
  return(res)
}


#' @importFrom cli cli_alert_warning
cds_warn <- function(cds) {
  if (!is.null(cds)) {
    cli::cli_alert_warning("The cds argument has been deprecated.  Passing this value to obj.")
    obj <- cds
  }
}

obj_stop <- function(obj) {
  stopifnot(
    "You must use this function on a cell_data_set or Seurat object" =
      class(obj) %in% c("cell_data_set", "Seurat", "monocle3", "SeuratObject")
  )
}
