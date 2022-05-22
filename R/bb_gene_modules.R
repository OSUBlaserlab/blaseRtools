#' @title A function to generate gene modules and add them to the Gene Metadata
#' @description Based on Monocle3's gene module functions.  Implemented with default values.  Will convert a Seurat object to a cell data set using SeuratWrappers and then calculate modules.  The function returns an object of the same type.
#' @param obj A single cell object of type Seurat or cell_data_set.
#' @param n_cores Number of processor cores to use for the analysis, Default: 8
#' @param cds Provided for backward compatibility.  If a value is supplied, it will return a warning and transfer to the obj argument., Default: NULL
#' @return An object of the same type:  Seurat or cell_data_set
#' @details see https://cole-trapnell-lab.github.io/monocle3/docs/differential/#gene-modules
#' @seealso
#'  \code{\link[SeuratWrappers]{as.cell_data_set}}
#'  \code{\link[monocle3]{graph_test}}, \code{\link[monocle3]{find_gene_modules}}
#'  \code{\link[dplyr]{rename}}, \code{\link[dplyr]{mutate-joins}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}
#'  \code{\link[forcats]{fct_shift}}
#' @rdname bb_gene_modules
#' @export
#' @importFrom SeuratWrappers as.cell_data_set
#' @importFrom monocle3 graph_test find_gene_modules
#' @importFrom dplyr rename left_join mutate select
#' @importFrom forcats fct_shift
bb_gene_modules <- function(obj,
                            n_cores = 8,
                            cds = NULL) {
  cds_warn(cds)
  obj_stop(obj)

  if ("Seurat" %in% class(obj)) {
    cds <- SeuratWrappers::as.cell_data_set(obj)
  } else {
    cds <- obj
  }

  graph_test_res <-
    monocle3::graph_test(
      cds = cds,
      neighbor_graph = "knn",
      cores = n_cores,
      verbose = TRUE
    )
  mod_deg_ids <- row.names(subset(graph_test_res, q_value < 0.05))

  gene_module_df <-
    monocle3::find_gene_modules(cds[mod_deg_ids,], cores = n_cores) |>
    dplyr::rename(feature_id = id)

  gene_module_df <-
    dplyr::left_join(bb_rowmeta(obj), gene_module_df, by = "feature_id") |>
    dplyr::mutate(module_labeled = ifelse(is.na(module), "No Module", paste0("Module ", module))) |>
    dplyr::mutate(module_labeled = factor(module_labeled, levels = str_sort(unique(module_labeled), numeric = TRUE))) |>
    dplyr::mutate(supermodule_labeled = ifelse(
      is.na(module),
      "No Supermodule",
      paste0("Supermodule ", module)
    )) |>
    dplyr::mutate(supermodule_labeled = factor(supermodule_labeled, levels = str_sort(unique(
      supermodule_labeled
    ), numeric = TRUE))) |>
    dplyr::mutate(supermodule_labeled = forcats::fct_shift(supermodule_labeled)) |>
    dplyr::select(feature_id,
                  module,
                  module_labeled,
                  supermodule,
                  supermodule_labeled)

  obj <- bb_tbl_to_rowdata(obj, min_tbl = gene_module_df)

  return(obj)
}
