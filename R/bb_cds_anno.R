#' @title Annotate a CDS by Label Transfer
#' @description Use this function to transfer cell labels from one single cell dataset to another.  If a cds is provided for either reference or query, it is converted to a seurat object and the labels are transferred to the query by anchor finding.  The assignments will be in the form of a new cds column with name predicted."column name from reference".  This should be unique on the first application of the function to a query dataset.  However, if running queries against more than 1 reference data set it is possible that you will unintentionally generate the same column name which would overwrite the first assignment column.  The function checks for this and aborts with the recommendation to supply a unique id to the unique_id parameter.
#' @param query_cds The single cell data set you wish to annotate.  Must be a CDS.
#' @param ref The reference single cell data set.  May be either a CDS or Seurat object.
#' @param transfer_col The column from the reference data set that provides the labels.
#' @param unique_id A unique identifier to add to the column with the transferred labels. Default is NULL but it is recommended to provide an informative label when annotating against more than one reference.  Default is NULL.
#' @return a CDS with two new cell metadata columns
#' @seealso
#'  \code{\link[cli]{cli_abort}}
#'  \code{\link[Seurat]{components}}, \code{\link[Seurat]{SCTransform}}, \code{\link[Seurat]{RunPCA}}, \code{\link[Seurat]{RunUMAP}}, \code{\link[Seurat]{FindTransferAnchors}}, \code{\link[Seurat]{MapQuery}}
#'  \code{\link[monocle3]{exprs}}
#'  \code{\link[SingleCellExperiment]{character(0)}}
#'  \code{\link[purrr]{reexports}}
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate-joins}}
#'  \code{\link[tibble]{as_tibble}}
#' @rdname bb_cds_anno
#' @export
#' @importFrom cli cli_abort
#' @importFrom Seurat CreateSeuratObject AddMetaData SCTransform RunPCA RunUMAP Idents FindTransferAnchors VariableFeatures MapQuery
#' @importFrom monocle3 exprs
#' @importFrom SingleCellExperiment colData
#' @importFrom purrr set_names
#' @importFrom dplyr select left_join
#' @importFrom tibble as_tibble
bb_cds_anno <- function (query_cds,
                         ref,
                         transfer_col,
                         unique_id = NULL) {
    if (!is.null(unique_id)) {
      unique_id <- paste0(unique_id, "_")
    }
    if (paste0("predicted.", unique_id, transfer_col) %in% colnames(bb_cellmeta(query_cds))) {
      cli::cli_abort(
        "This function will overwrite existing data in your cds. Provide a new value for unique_id.  Aborting."
      )
    }
    if ("cell_data_set" %in% class(ref)) {
      ref_seurat <- Seurat::CreateSeuratObject(counts = monocle3::exprs(ref))
      ref_seurat <-
        Seurat::AddMetaData(ref_seurat,
                    SingleCellExperiment::colData(ref)[[transfer_col]],
                    col.name = transfer_col) |>
        Seurat::SCTransform() |>
        Seurat::RunPCA() |>
        Seurat::RunUMAP(dims = 1:30)
    }
    else if ("Seurat" %in% class(ref)) {
      ref_seurat <- ref
    }
    else {
      return("The reference must be either a cell data set or Seurat object.")
    }
    Seurat::Idents(ref_seurat) <- transfer_col
    if ("cell_data_set" %in% class(query_cds)) {
      query_seurat <- Seurat::CreateSeuratObject(counts = monocle3::exprs(query_cds)) |>
        Seurat::SCTransform() |>
        Seurat::RunPCA()

    }
    else {
      return("The query must be a cell data set object.")
    }
    transfer_anchors <- Seurat::FindTransferAnchors(
      reference = ref_seurat,
      query = query_seurat,
      features = Seurat::VariableFeatures(object = ref_seurat),
      reference.reduction = "pca"
    )
    query_ref <-
      Seurat::MapQuery(
        anchorset = transfer_anchors,
        query = query_seurat,
        reference = ref_seurat,
        refdata = list(transfer_col = transfer_col)
      )
    # return(query_ref)
    cols_to_add <- purrr::set_names(
      dplyr::select(
        dplyr::left_join(
          bb_cellmeta(query_cds),
          tibble::as_tibble(query_ref@meta.data, rownames = "cell_id"),
          by = "cell_id"
        ),
        cell_id,
        predicted.transfer_col,
        predicted.transfer_col.score
      ),
      c(
        "cell_id",
        paste0("predicted.", unique_id, transfer_col),
        paste0("predicted.", unique_id, transfer_col, ".score")
      )
    )
    result <-
      bb_tbl_to_coldata(obj = query_cds, min_tbl = cols_to_add)
    return(result)
  }

