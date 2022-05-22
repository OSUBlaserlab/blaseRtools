#' A Function To Add Tibble Columns To Cell Metadata
#'
#' @param obj A Seurat or cell data set object
#' @param min_tbl A tibble containing only the columns you want to add plus one column for joining.  Cell IDs may not be duplicated but missing cells are ok; values will be replaced by NA.
#' @param join_col The column in min_tbl containing the join information for the cds rowData. Defaults to "cell_id".
#' @param cds Retained for backwards compatibility.  If supplied, will generate a warning and pass argument to obj.  Default = NULL
#' @return An object of the same class
#' @export
#' @rdname bb_tbl_to_coldata
#' @export
#' @importFrom dplyr select rename left_join
#' @importFrom tidyselect contains
#' @importFrom SummarizedExperiment colData
bb_tbl_to_coldata <-
  function (obj,
            min_tbl,
            join_col = "cell_id",
            cds = NULL) {
    cds_warn(cds)
    obj_stop(obj)
    cell_ids <- bb_cellmeta(obj) |>
      dplyr::select(cell_id)
    min_tbl <-
      dplyr::rename(min_tbl, cell_id = tidyselect::contains(join_col))
    full_tbl <-
      dplyr::left_join(cell_ids,
                       min_tbl,
                       by = "cell_id")
    stopifnot(
      "Joining min_tbl to colData failed.  Check for duplicates in min_tbl." =
        identical(full_tbl$cell_id, bb_cellmeta(obj)$cell_id)
    )
    cols_to_add <- full_tbl |>
      dplyr::select(-cell_id)

    if ("cell_data_set" %in% class(obj)) {
      for (i in 1:ncol(cols_to_add)) {
        SummarizedExperiment::colData(obj)$new <- cols_to_add |> pull(i)
        names(SummarizedExperiment::colData(obj))[names(SummarizedExperiment::colData(obj)) == "new"] <-
          colnames(cols_to_add)[i]

      }
    } else {
      for (i in 1:ncol(cols_to_add)) {
        obj@meta.data[colnames(cols_to_add)[i]] <- cols_to_add |> pull(i)
      }
    }
    return(obj)
  }
