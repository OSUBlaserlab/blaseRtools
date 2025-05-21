#' Rejoin Blinded Scores to Original File Names
#'
#' @description Will rejoin scoresheet and blinded key to produce unblinded results.  If you change the names of either of those files, they have to be provided as arguments to the function.  Otherwise keyfile and scorefile are optional.
#' @param directory The linux-style filepath of the folder containing the scoresheet and blinded key.
#' @param keyfile Optional:  filename of the key file.  Defaults to "blinding_key.csv".
#' @param scorefile Optional:  filename of the score file.  Defaults to "scoresheet.csv".
#' @param analysis_file Complete file path to the the unblinded main analysis sheet.  The function will will left_join analysis_file and unblinded results.  In the process, it will necessarily convert windows file paths to linux-style file paths.  Samples not included in the blinding should return with NA values for the added columns.  New data columns being added on from scoresheet should be unique relative to analysis_file.
#' @param file_column The column in analysis_file containing file paths for the files that were blinded.
#' @return nothing
#' @export
#' @import tidyverse fs
bb_unblind_images <-
  function(directory,
           keyfile = "blinding_key.csv",
           scorefile = "scoresheet.csv",
           analysis_file,
           file_column)  {
    key <- readr::read_csv(fs::path(directory, keyfile))
    scores <- readr::read_csv(fs::path(directory, scorefile))
    res <- dplyr::left_join(key,
                     scores) |>
      dplyr::arrange(source_file)
    readr::write_csv(res, file = fs::path(directory, "unblinded_result.csv"))
    analysis_file <-
      readr::read_csv(analysis_file) |>
      dplyr::mutate(dplyr::across(tidyr::contains(file_column), bb_fix_file_path)) |>
      dplyr::rename(source_file = tidyr::contains(file_column))
    analysis_joined <-
      dplyr::left_join(analysis_file, res) |>
      dplyr::rename(!!file_column := source_file) |>
      dplyr::select(-blinded_file)
    return(analysis_joined)

  }
