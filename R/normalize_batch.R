#' @title Noramalize a Data Table by Group and Batch
#' @description Often you will have a data table with repeats or batches of the same experiment.  An effective way to control for batch effects is to normalize the data from each batch to a control group present in all of the experiments.  To use this function, provide such a data table, identify the column holding the experimental group data, the identity of the control group to normalize by, the column holding the batch data, and the column holding the numerical data to normalize.  Also select the function to average by (mean or median).  The function will return the data table with three new columns:  the average of the control group by batch, fold change of each observation relative to the batch average and the log2-transformed fold change.  This function is pipe-friendly.
#' @param data a tibble
#' @param group_col the column containing the experimental group identifier
#' @param norm_group the experimental group you want to normalize to across batches
#' @param batch_col the column containing the batch identifier
#' @param data_col the column with your data
#' @param fun averaging function to use, Default: c("mean", "median")
#' @return A tibble with new columns indicating batch normalization group average, fold change for each observation relative to the batch average and log2 fold change
#' @seealso
#'  \code{\link[rlang]{arg_match}}
#'  \code{\link[cli]{cli_abort}}
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{group_by}}, \code{\link[dplyr]{summarise}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate-joins}}, \code{\link[dplyr]{mutate}}
#' @rdname normalize_batch
#' @export
#' @importFrom rlang arg_match
#' @importFrom cli cli_abort
#' @importFrom dplyr filter group_by summarize select left_join mutate
normalize_batch <- function(data,
                            group_col,
                            norm_group,
                            batch_col,
                            data_col,
                            fun = c("mean", "median")) {
  # check the arguments
  fun <- rlang::arg_match(fun)
  if (class(data[[group_col]]) %notin% c("character", "factor", "integer"))
    cli::cli_abort("The group_col argument must be either a character, factor or integer column.")
  if (class(data[[batch_col]]) %notin% c("character", "factor", "integer"))
    cli::cli_abort("The batch_col argument must be either a character, factor or integer column.")
  if (norm_group %notin% data[[group_col]])
    cli::cli_abort("The normalization group is not present in group_col.")


  # calculate the normalization values
  norm_vals <- data |>
    dplyr::filter(!!sym(group_col) == norm_group) |>
    dplyr::group_by(!!sym(batch_col)) |>
    dplyr::summarize(batch_mean = mean(!!sym(data_col)),
              batch_median = median(!!sym(data_col))) |>
    dplyr::select(!!sym(batch_col), !!sym(paste0("batch_", fun)))

  # add the normalization values back onto the data by group

  data <- data |>
    dplyr::left_join(norm_vals, by = batch_col) |>
    dplyr::mutate(fold_change_batch = !!sym(data_col)/!!sym(paste0("batch_", fun))) |>
    dplyr::mutate(l2fc_batch = log2(fold_change_batch))
  data
}
