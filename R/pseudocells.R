#' @title Generate a CDS from Bulk RNA-seq data
#' @description This function allows you to simulate single cell data from bulk RNA-seq data.  It requires a TPM count matrix.  The rationale is that TPM quantifies transcript counts per million reads, so you can think of this like 1 million UMI counts from a scRNA-seq experiment distributed in a certain way across the transcriptome.  This function samples *n_pseudocells*  with *transcripts_per_pseudocell* from this distribution.  Then it creates a cell data set based on a matrix of these samples.
#'
#' Importantly, this function cannot accurately identify transcriptional heterogenity within the bulk data.  The sampling effect may reveal some potential heterogeneity but there is no method here for determining whether this is due to randomness or heterogeneity within the data.
#'
#' The intended use of this function is to project a pseudocell cds onto an actual single cell data set.  This can be used to help identify regions of UMAP space that are shared between the pseudocells and the real cells.
#'
#' Note:  The matrix rownames and row metadata for the returned cds will contain only the rownames from the input tpm_matrix.  So these should share the same namespace as the single cell cds that the pseudocell data will be projected onto.
#'
#' @param tpm_matrix A matrix of TPM counts from a bulk RNA experiment.  A cds will be generated for each column in this matrix and the result combined.
#' @param n_pseudocells The number of pseudocells to create.  Should be length 1 or to specify a uniqe n_pseudocells for each dataset, a vector of the same length as the number of columns in tpm_matrix.  This value will be recycled if necessary to match the number of columns in tpm_mtx.
#' @param transcripts_per_pseudocell The number of transcripts to sample for each pseudocell.  Should be similar to the median number of UMI in the single cell data the pseudocells will be projected onto.  Should be length 1 or to specify a uniqe n_pseudocells for each dataset, a vector of the same length as the number of columns in tpm_matrix.  This value will be recycled if necessary to match the number of columns in tpm_mtx.
#' @param remove_genes A vector of genes to remove before sampling.  Should be the same or similar to the genes removed from the single cell data the pseudocells will be projected onto, Default: NULL
#' @return A cell data set
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[purrr]{map}}, \code{\link[purrr]{reduce}}, \code{\link[purrr]{map2}}, \code{\link[purrr]{pmap}}
#'  \code{\link[tibble]{tibble}}
#'  \code{\link[dplyr]{count}}, \code{\link[dplyr]{select}}
#'  \code{\link[Seurat]{components}}
#'  \code{\link[monocle3]{new_cell_data_set}}, \code{\link[monocle3]{combine_cds}}, \code{\link[monocle3]{preprocess_cds}}, \code{\link[monocle3]{reduce_dimension}}
#' @rdname bb_pseudocells
#' @export
#' @importFrom cli cli_abort cli_warn
#' @importFrom purrr pmap map reduce map2
#' @importFrom tibble tibble
#' @importFrom dplyr count select
#' @importFrom Seurat as.sparse
#' @importFrom monocle3 new_cell_data_set combine_cds preprocess_cds reduce_dimension
bb_pseudocells <-
  function(tpm_matrix,
           n_pseudocells,
           transcripts_per_pseudocell,
           remove_genes = NULL) {
    if (!is.numeric(n_pseudocells)) {
      cli::cli_abort("{.emph n_pseudocells} must be a numeric.")
    }
    if (!is.numeric(transcripts_per_pseudocell)) {
      cli::cli_abort("{.emph transcripts_per_pseudocell} must be a numeric.")
    }

    # if (ncol(tpm_matrix) != length(n_pseudocells)) {
    #   cli::cli_warn(
    #     "Recycling {.emph n_pseudocells} to match the number of columns in {.emph tpm_matrix}"
    #   )
    # }
    #
    # if (ncol(tpm_matrix) != length(transcripts_per_pseudocell)) {
    #   cli::cli_warn(
    #     "Recycling {.emph transcripts_per_pseudocell} to match the number of columns in {.emph tpm_matrix}"
    #   )
    # }
    #

    exprs_list <-
      purrr::pmap(
        .l = list(
          sample_name = colnames(tpm_matrix),
          np = n_pseudocells,
          tpp = transcripts_per_pseudocell
        ),
        .f = \(sample_name,
               np,
               tpp,
               dat = tpm_matrix,
               remove = remove_genes) {
          bulk <- sample_name
          vec <- dat[, bulk]
          vec <- vec[!names(vec) %in% remove]
          vec <- rep(names(vec), times = vec)

          count_list <- purrr::map(.x = 1:np,
                                   .f = \(
                                     x,
                                     sample_size = tpp,
                                     cell = bulk,
                                     v = vec
                                   ) {
                                     cell <- paste0(cell, "_pc_", x)
                                     v <-
                                       sample(v,
                                              size = sample_size,
                                              replace = FALSE)
                                     tibble::tibble(original_rowname = v) |>
                                       dplyr::count(original_rowname, name = cell)
                                   })

          full_tbl <- purrr::reduce(count_list,
                                    .f = full_join,
                                    by = "original_rowname")
          full_tbl[is.na(full_tbl)] <- 0

          full_mtx <- as.matrix(full_tbl |>
                                  dplyr::select(-original_rowname))
          rownames(full_mtx) <-
            full_tbl$original_rowname
          Seurat::as.sparse(full_mtx)

        }
      ) |> set_names(colnames(tpm_matrix))


    coldata_list <- purrr::map2(.x = exprs_list,
                                .y = colnames(tpm_matrix),
                                .f = \(x, y) {
                                  coldata <- tibble::tibble(cell_id = colnames(x),
                                                            dataset = y) |>
                                    base::as.data.frame()
                                  rownames(coldata) <- coldata$cell_id
                                  coldata$cell_id <- NULL
                                  coldata
                                })

    rowdata_list <- purrr::map(.x = exprs_list,
                               .f = \(x) {
                                 data.frame(
                                   id = rownames(x),
                                   gene_short_name = rownames(x),
                                   row.names = rownames(x)
                                 )

                               })



    cell_line_cds_list <- purrr::pmap(
      .l = list(expr = exprs_list,
                col = coldata_list,
                row = rowdata_list),
      .f = \(expr, col, row) {
        monocle3::new_cell_data_set(
          expression_data = expr,
          cell_metadata = col,
          gene_metadata = row
        )
      }
    )


    cell_line_cds <-
      monocle3::combine_cds(cds_list = cell_line_cds_list)
    cell_line_cds <- monocle3::preprocess_cds(cell_line_cds)
    cell_line_cds <- monocle3::reduce_dimension(cell_line_cds)
    cell_line_cds

  }
