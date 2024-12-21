#' @title Annotate Single cell Data Using Monocle3
#' @description This function wraps the monocle3 method for data projection and label transfer into a function.  This function takes in reference and query cds objects plus a character vector of label column names and returns the query CDS with the new labels.
#' @param cds_qry A query cell data set.
#' @param cds_ref A reference cell data set.
#' @param labels A character vector of cell metadata column names to transfer.
#' @param suffix A character string of length 1 to append to all of the tranferred column names.  There is no checking for name conflicts, so use this sensibly to prevent overwriting preexisting columns.  Default = "_ref".
#' @return A cell data set.
#' @details see https://cole-trapnell-lab.github.io/monocle3/docs/projection/
#' @seealso
#'  \code{\link[monocle3]{estimate_size_factors}}, \code{\link[monocle3]{preprocess_cds}}, \code{\link[monocle3]{reduce_dimension}}, \code{\link[monocle3]{save_transform_models}}, \code{\link[monocle3]{load_transform_models}}, \code{\link[monocle3]{preprocess_transform}}, \code{\link[monocle3]{reduce_dimension_transform}}, \code{\link[monocle3]{transfer_cell_labels}}, \code{\link[monocle3]{fix_missing_cell_labels}}
#'  \code{\link[fs]{path}}
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{reexports}}
#' @rdname bb_monocle_anno
#' @export
#' @importFrom monocle3 estimate_size_factors preprocess_cds reduce_dimension save_transform_models load_transform_models preprocess_transform reduce_dimension_transform transfer_cell_labels fix_missing_cell_labels
#' @importFrom fs path
#' @importFrom dplyr select matches
#' @import cli
#' @importFrom function purrr reduce map
#' @importFrom dplyr select join_by left_join
bb_monocle_anno <-
  function(cds_qry, cds_ref, labels, suffix = "_ref") {
    cds_stop(cds_qry)
    cds_stop(cds_ref)
    cds_qry_hold <- cds_qry
    cols <- colnames(bb_cellmeta(cds_ref))

    if (!all(labels %in% cols)) {
      cli::cli_abort("All labels must be column names in cds_ref.")
    }

    # Genes in reference.
    genes_ref <- row.names(cds_ref)

    # Genes in query.
    genes_qry <- row.names(cds_qry)

    # Shared genes.
    genes_shared <- intersect(genes_ref, genes_qry)

    if (!all(genes_qry == genes_shared)) {
      cli::cli_alert_info("Reprocessing cds_qry...")
      cds_qry <- cds_qry[genes_shared,]
      cds_qry <- monocle3::estimate_size_factors(cds_qry)
      cli::cli_alert_info("...Done")
    }

    if (!all(genes_ref == genes_shared)) {
      cli::cli_alert_info("Reprocessing cds_ref...")
      cds_ref <- cds_ref[genes_shared,]
      cds_ref <- monocle3::estimate_size_factors(cds_ref)
      cds_ref <- monocle3::preprocess_cds(cds_ref, num_dim = 100)
      cds_ref <-
        monocle3::reduce_dimension(cds_ref, build_nn_index = TRUE)
      cli::cli_alert_info("...Done")
    }
    # Save the PCA and UMAP transform models for use with projection.

    monocle3::save_transform_models(cds_ref, fs::path(tempdir(), "cds_ref_test_models"))

    # Load the reference transform models into the query cds.
    cli::cli_alert_info("Loading cds_ref models and reprocessing cds_qry...")
    cds_qry <-
      monocle3::load_transform_models(cds_qry, fs::path(tempdir(), "cds_ref_test_models"))

    # Apply the reference transform models to the query cds.
    cds_qry <- monocle3::preprocess_transform(cds_qry)
    cds_qry <- monocle3::reduce_dimension_transform(cds_qry)
    cli::cli_alert_info("...Done")

    new_cellmeta_list <- purrr::map(.x = labels,
                             .f = \(x,
                                    qry = cds_qry,
                                    ref = cds_ref,
                                    suf = suffix) {
                               cds_qry_lab_xfr <-
                                 monocle3::transfer_cell_labels(
                                   qry,
                                   reduction_method = "UMAP",
                                   ref_coldata = colData(ref),
                                   ref_column_name = x,
                                   query_column_name = paste0(x, "_xfr"),
                                   transform_models_dir = fs::path(tempdir(), "cds_ref_test_models")
                                 )

                               cds_qry_lab_fix <-
                                 monocle3::fix_missing_cell_labels(
                                   cds_qry_lab_xfr,
                                   reduction_method = "UMAP",
                                   from_column_name = paste0(x, "_xfr"),
                                   to_column_name = paste0(x, suf)
                                 )
                               bb_cellmeta(cds_qry_lab_fix) |>
                                 dplyr::select(cell_id, paste0(x, suf))

                             })
    cli::cli_alert_success("Transferred labels to cds_qry.")
    new_cellmeta <-
      purrr::reduce(.x  = new_cellmeta_list,
                    .f  = dplyr::left_join,
                    by = dplyr::join_by("cell_id"))

    bb_tbl_to_coldata(cds_qry_hold, min_tbl = new_cellmeta)



  }
