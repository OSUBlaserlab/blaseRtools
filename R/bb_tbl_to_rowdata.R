#' A Function To Add Tibble Columns To Feature Metadata
#'
#' @param obj A Seurat or cell data set object
#' @param assay The assay to which to add the metadata column, Default = RNA
#' @param min_tbl A tibble containing only the columns you want to add plus one column for joining.  Features cannot be duplicated but missing features are ok and will be replaced by NA.
#' @param join_col The column in min_tbl containing the join information for the cds rowData. Defaults to "feature_id".
#' @param cds Retained for backwards compatibility.  If supplied, will generate a warning and pass argument to obj.  Default = NULL
#' @return An object of the same class
#' @rdname bb_tbl_to_rowdata
#' @export
#' @importFrom dplyr select rename left_join
#' @importFrom tidyselect contains
#' @importFrom SingleCellExperiment rowData
#' @importFrom SeuratObject AddMetaData
bb_tbl_to_rowdata <-
  function (obj,
            assay = "RNA",
            min_tbl,
            join_col = "feature_id",
            cds = NULL) {
    blaseRtools:::cds_warn(cds)
    blaseRtools:::obj_stop(obj)
    row_ids <- bb_rowmeta(obj) |>
      dplyr::select(feature_id)
    min_tbl <-
      dplyr::rename(min_tbl, feature_id = tidyselect::contains(join_col))
    full_tbl <-
      dplyr::left_join(row_ids,
                       min_tbl,
                       by = "feature_id")
    stopifnot(
      "Joining min_tbl to rowData failed.  Check for duplicates in min_tbl." =
        identical(full_tbl$feature_id, bb_rowmeta(obj)$feature_id)
    )
    cols_to_add <- full_tbl |>
      dplyr::select(-feature_id)

    if ("cell_data_set" %in% class(obj)) {
      for (i in 1:ncol(cols_to_add)) {
        SummarizedExperiment::rowData(obj)$new <- cols_to_add |> pull(i)
        names(SummarizedExperiment::rowData(obj))[names(SummarizedExperiment::rowData(obj)) == "new"] <-
          colnames(cols_to_add)[i]

      }
    } else {
      for (i in 1:ncol(cols_to_add)) {
        obj[[assay]] <- SeuratObject::AddMetaData(object = obj[[assay]],
                                                  metadata = cols_to_add[i])
      }


    }
    return(obj)
  }
