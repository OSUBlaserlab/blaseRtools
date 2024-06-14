#' @title Create a Loupe File from a Cell Data Set
#' @description Converts a cell data set into a loupe file.
#' @param cds Input cds.  Only works on cell data set objects.  For seurat objects, use built in loupeR functions
#' @param output_dir Output directory, Default: '.'
#' @param output_file Name of the loupe file.  .cloupe will be appended for compatibility, Default: 'loupe'
#' @return Nothing
#' @seealso
#'  \code{\link[cli]{cli_abort}}, \code{\link[cli]{cli_alert}}
#'  \code{\link[fs]{create}}
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{across}}, \code{\link[dplyr]{reexports}}
#'  \code{\link[tidyselect]{starts_with}}
#'  \code{\link[SingleCellExperiment]{reducedDims}}
#'  \code{\link[monocle3]{exprs}}
#'  \code{\link[waldo]{compare}}
#'  \code{\link[loupeR]{create_loupe}}
#' @rdname bb_loupeR
#' @export
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom fs dir_create
#' @importFrom dplyr select mutate across everything
#' @importFrom tidyselect matches
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom monocle3 exprs
#' @importFrom waldo compare
#' @importFrom loupeR create_loupe
bb_loupeR <- function(cds,
                      output_dir = ".",
                      output_file = "loupe") {
  if (class(cds) != "cell_data_set") {
    cli::cli_abort("This function only works with a cell_data_set.\nFor Seurat functions, use loupeR::create_loupe_from_seurat.")
  }

  fs::dir_create(output_dir)

  cli::cli_alert_info("Dropping numeric and logical cell metadata variables.\n")

  cellmeta <- bb_cellmeta(cds) |>
    dplyr::select(!where(is.numeric)) |>
    dplyr::select(!where(is.logical)) |>
    dplyr::select(-tidyselect::matches("cell_id", "barcode")) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.factor)) |>
    as.list()

  proj <- list(SingleCellExperiment::reducedDims(cds)[["UMAP"]])
  names(proj) <- "umap"

  mat <- monocle3::exprs(cds)
  waldo::compare(rownames(proj$umap), colnames(mat))

  cli::cli_alert_info("Creating loupe file:\n")
  loupeR::create_loupe(
    count_mat = mat,
    clusters = cellmeta,
    projections = proj,
    output_dir = output_dir,
    output_name = output_file
  )

}
