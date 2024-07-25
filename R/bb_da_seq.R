#' @title Calculate Differential Abundance Using DAseq
#' @description This function takes a cell data set object and a cell metadata variable as input.  The latter is specified as a named list.  The name is the cell metadata variable name and each element must be a character vector of length 2 specifying the levels of that variable to compare.  For example, for a column named "genotype" with levels "WT", "heterozygote", and "homozygote", you would compare differential abundance between WT and homozygote like this:  list(genotype = c("WT", "homozygote")).  Cells with these classifications will get a differential abundance score.  Negative values will be assigned to the first element (WT in this example) and positive values to the second element (homozygote in this example).  Cells with other values for this variable ("heterozygote" in this example) will get an NA.  Multiple comparisons can be performed for different levels within the same variable or for different variables altogether by providing additional elements to the list.  E.g. list(genotype = c("WT", "heterozygote), genotype = c("WT", "homozygote")).  For each comparison, a new cell metadata column will be added in the form of da_score_name_level1_level2.  E.g. da_score_genotype_WT_homozygote.
#' @param obj a cell data set object
#' @param comparison_list a named list as specified in description
#' @param sample_var the cell metadata variable holding biological sample information, Default: 'sample'
#' @return a cell data set
#' @seealso
#'  \code{\link[cli]{cli_div}}, \code{\link[cli]{cli_alert}}
#'  \code{\link[purrr]{map2}}, \code{\link[purrr]{reexports}}, \code{\link[purrr]{pmap}}, \code{\link[purrr]{reduce}}
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{pull}}, \code{\link[dplyr]{mutate-joins}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{join_by}}, \code{\link[dplyr]{rename}}
#'  \code{\link[SingleCellExperiment]{reducedDims}}
#'  \code{\link[DAseq]{getDAcells}}
#'  \code{\link[tibble]{as_tibble}}
#' @rdname bb_daseq
#' @export
#' @importFrom cli cli_div cli_alert_info
#' @importFrom purrr map2 set_names pmap reduce
#' @importFrom dplyr filter pull left_join select join_by rename_with
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom DAseq getDAcells
#' @importFrom tibble as_tibble
bb_daseq <- function(obj, comparison_list, sample_var = "sample") {
  obj_stop(obj)
  cli::cli_div(theme = list(span.emph = list(color = "orange")))
  da_cells <- purrr::map2(.x = comparison_list,
                   .y = names(comparison_list),
                   .f = \(x,
                          y,
                          dat = obj,
                          svar = sample_var) {
                     cli::cli_alert_info("Calculating differential abundance based on the {.emph {y}} variable.")
                     cli::cli_alert_info(
                       "Cells from the {.emph {x[1]}} group will be assigned negative values and cells from the {.emph {x[2]}} group will be assigned positive values."
                     )
                     cli::cli_alert_info("Cells from other groups in this variable will return {.emph {NA}} with a warning.")
                     labels_1 <- bb_cellmeta(dat) |>
                       dplyr::filter(.data[[y]] == x[1]) |>
                       dplyr::pull(.data[[svar]]) |>
                       unique()
                     labels_2 <- bb_cellmeta(dat) |>
                       dplyr::filter(.data[[y]] == x[2]) |>
                       dplyr::pull(.data[[svar]]) |>
                       unique()
                     pcas <-
                       SingleCellExperiment::reducedDims(obj)$PCA[, 1:30]
                     umaps <-
                       SingleCellExperiment::reducedDims(obj)$UMAP
                     cell_labels <-
                       bb_cellmeta(obj) |>
                       dplyr::pull(.data[[svar]])
                     DAseq::getDAcells(
                       X = pcas,
                       cell.labels = cell_labels,
                       labels.1 = labels_1,
                       labels.2 = labels_2,
                       k.vector = seq(50, 500, 50),
                       plot.embedding = umaps
                     )



                   }) |>
    purrr::set_names(names(comparison_list))

  tbl_to_add <- purrr::pmap(.l = list(
    x = da_cells,
    y = names(comparison_list),
    z = comparison_list
  ),
  .f = \(x,
         y,
         z,
         dat = obj) {
    new_column <- paste0("da_score_", y, "_", z[1], "_", z[2])
    cli::cli_alert_info("Adding differential abundance to a new cell metadata column, {.emph {new_column}}")
    dplyr::left_join(
      bb_cellmeta(dat),
      tibble::as_tibble(x$pred.plot[["data"]], rownames = "cell_id") |>
        dplyr::select(cell_id, da_score = Score),
      by = dplyr::join_by("cell_id")
    ) |>
      dplyr::select(cell_id, da_score) |>
      dplyr::rename_with(.fn = ~ paste0("da_score_", y, "_", z[1], "_", z[2]),
                  .cols = da_score)

  }) |>
    purrr::reduce(.f = dplyr::left_join, by = dplyr::join_by(cell_id))

  bb_tbl_to_coldata(obj, min_tbl = tbl_to_add)


}
