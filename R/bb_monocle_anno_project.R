#' @title Annotate Single cell Data Using Monocle3
#' @description This function wraps the monocle3 method for data projection and label transfer into a function.  This function takes in reference and query cds objects plus a character vector of label column names and returns the query CDS with the new labels.


#' @title Monocle Annotation and Projections
#' @description bb_monocle_anno and bb_monocle_project use similar inputs and methods and both return cds objects, but the cds objects are different.  Both wrap around the monocle3 mehtod for data projection and label transfer.
#'
#' bb_monocle_anno takes in reference and query cds objects plus a character vector of label column names to transfer.  It returns the query CDS with new cell metadata contining one or more columns with the labels corresponding to the nearest neighbor in the reference data.  A suffix is appended to the new column name in the result. By default this is "_ref" but it can be changed using the suffix parameter.
#'
#' bb_monocle_project also takes in reference and query cds objects and character vector of label column names.  In this case, a combined cds object is returned carrying the query and reference data projected int the reference data space.  New column names are added indicating the reference and query data.  The query cells are given the label of their nearest neighbor in the reference.  These label column(s) are prepended with "merged_".
#'
#'Importantly, for both of these to work, the reference and query genes must have a shared namespace.
#'
#'
#' @param cds_qry A query cell data set.
#' @param cds_ref A reference cell data set.
#' @param labels A character vector of cell metadata column names to transfer.
#' @param suffix A character string of length 1 to append to all of the tranferred column names.  There is no checking for name conflicts, so use this sensibly to prevent overwriting preexisting columns.  Default = "_ref".
#' @param use_aligned Whether to use aligned PCA coordinates from cds_ref, if they are available.  Default = TRUE.
#' @return A cell data set
#' @details see https://cole-trapnell-lab.github.io/monocle3/docs/projection/
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[cli]{cli_abort}}, \code{\link[cli]{cli_alert}}
#'  \code{\link[base]{sets}}
#'  \code{\link[monocle3]{estimate_size_factors}}, \code{\link[monocle3]{preprocess_cds}}, \code{\link[monocle3]{reduce_dimension}}, \code{\link[monocle3]{save_transform_models}}, \code{\link[monocle3]{load_transform_models}}, \code{\link[monocle3]{preprocess_transform}}, \code{\link[monocle3]{reduce_dimension_transform}}, \code{\link[monocle3]{transfer_cell_labels}}, \code{\link[monocle3]{fix_missing_cell_labels}}
#'  \code{\link[SingleCellExperiment]{reducedDims}}
#'  \code{\link[fs]{path}}
#'  \code{\link[purrr]{map}}, \code{\link[purrr]{reduce}}
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate-joins}}, \code{\link[dplyr]{join_by}}
#' @rdname bb_monocle_
#' @importFrom cli cli_abort cli_alert_info cli_alert_success
#' @importFrom base intersect
#' @importFrom monocle3 estimate_size_factors preprocess_cds reduce_dimension save_transform_models load_transform_models preprocess_transform reduce_dimension_transform transfer_cell_labels fix_missing_cell_labels
#' @importFrom SingleCellExperiment reducedDimNames reducedDim reducedDims
#' @importFrom fs path
#' @importFrom purrr map reduce
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr select left_join join_by
bb_monocle_ <-
  function(cds_qry,
           cds_ref,
           labels,
           suffix,
           use_aligned) {
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
    genes_shared <- base::intersect(genes_ref, genes_qry)

    if (!all(genes_qry == genes_shared)) {
      cli::cli_alert_info("Reprocessing cds_qry...")
      cds_qry <- cds_qry[genes_shared,]
      cds_qry <- monocle3::estimate_size_factors(cds_qry)
      cli::cli_alert_info("...Done")
    }

    cli::cli_alert_info("Reprocessing cds_ref...")
    if ("Aligned" %in% SingleCellExperiment::reducedDimNames(cds_ref)) {
      if (use_aligned) {
        cli::cli_alert_info("Using aligned PCA values for the projection...")
        # replace the PCA values
        SingleCellExperiment::reducedDim(cds_ref) <-
          SingleCellExperiment::reducedDims(cds_ref)$Aligned

      } else {
        cli::cli_alert_info("Ignoring alignment and using original PCA values...")
      }
    }
    cds_ref <- cds_ref[genes_shared, ]
    cds_ref <- monocle3::estimate_size_factors(cds_ref)
    cds_ref <- monocle3::preprocess_cds(cds_ref, num_dim = 100)
    cds_ref <-
      monocle3::reduce_dimension(cds_ref,
                                 preprocess_method = "PCA",
                                 build_nn_index = TRUE)
    cli::cli_alert_info("...Done")
    # Save the PCA and UMAP transform models for use with projection.

    monocle3::save_transform_models(cds_ref,
                                    fs::path(tempdir(), "cds_ref_test_models"),
                                    verbose = FALSE)

    # Load the reference transform models into the query cds.
    cli::cli_alert_info("Loading cds_ref models and reprocessing cds_qry...")
    cds_qry <-
      monocle3::load_transform_models(cds_qry, fs::path(tempdir(), "cds_ref_test_models"))

    # Apply the reference transform models to the query cds.
    cds_qry <- monocle3::preprocess_transform(cds_qry)
    cds_qry <- monocle3::reduce_dimension_transform(cds_qry, )
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
                                          ref_coldata = SummarizedExperiment::colData(ref),
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
      purrr::reduce(
        .x  = new_cellmeta_list,
        .f  = dplyr::left_join,
        by = dplyr::join_by("cell_id")
      )

    cds_qry_final <-
      bb_tbl_to_coldata(cds_qry_hold, min_tbl = new_cellmeta)



    return(
      list(
        cds_qry_final = cds_qry_final,
        cds_ref = cds_ref,
        cds_qry = cds_qry,
        new_cellmeta = new_cellmeta
      )
    )
  }

#' @rdname bb_monocle_
#' @export
bb_monocle_anno <-
  function(cds_qry,
           cds_ref,
           labels,
           suffix = "_ref",
           use_aligned = TRUE) {
    monocle_return <- bb_monocle_(
      cds_qry = cds_qry,
      cds_ref = cds_ref,
      labels = labels,
      suffix = suffix,
      use_aligned = use_aligned
    )
    monocle_return$cds_qry_final
  }

#' @rdname bb_monocle_
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom monocle3 combine_cds
#' @importFrom purrr map reduce
#' @importFrom dplyr select sym all_of bind_rows
bb_monocle_project <-
  function(cds_qry,
           cds_ref,
           labels,
           suffix = "_ref",
           use_aligned = TRUE) {
    monocle_return <- bb_monocle_(
      cds_qry = cds_qry,
      cds_ref = cds_ref,
      labels = labels,
      suffix = suffix,
      use_aligned = use_aligned
    )
    # combine
    cds_ref <- monocle_return$cds_ref
    cds_qry <- monocle_return$cds_qry


    new_colname <- "data_set"
    disambiguator <- 1
    while (new_colname %in% c(colnames(colData(cds_ref)),
                              colnames(colData(cds_qry)))) {
      new_colname <-
        paste0(str_remove(new_colname, "\\..*"), ".", disambiguator)
      disambiguator <- disambiguator + 1
    }

    SummarizedExperiment::colData(cds_ref)[[new_colname]] <-
      'reference'
    SummarizedExperiment::colData(cds_qry)[[new_colname]] <- 'query'

    # Combine the reference and query data sets.
    cds_combined <-
      monocle3::combine_cds(
        list(cds_ref, cds_qry),
        keep_all_genes = TRUE,
        cell_names_unique = TRUE,
        keep_reduced_dims = TRUE
      )


    new_cellmeta <- monocle_return$new_cellmeta

    # merge the transferred colData columns
    new_cellmeta <- purrr::map(labels, .f = \(x,
                                              dat = new_cellmeta,
                                              suf = suffix) {
      dat |>
        dplyr::select(cell_id,
                      !!dplyr::sym(paste0("merged_", x)) := dplyr::all_of(paste0(x, suffix)))
    }) |>
      purrr::reduce(.f = left_join)

    ref_cellmeta <- purrr::map(labels, .f = \(x,
                                              dat = bb_cellmeta(cds_ref),
                                              suf = suffix) {
      dat |>
        dplyr::select(cell_id,
                      !!dplyr::sym(paste0("merged_", x)) := dplyr::all_of(x))
    }) |>
      purrr::reduce(.f = left_join)
    new_cellmeta <- dplyr::bind_rows(new_cellmeta, ref_cellmeta)
    cds_combined <-
      bb_tbl_to_coldata(cds_combined, min_tbl = new_cellmeta)


  }

