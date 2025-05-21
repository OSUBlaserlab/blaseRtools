#' Make a Copy of Image Files and Rename With File Hashes in Blinded Folder
#'
#' @description Will copy and rename  the files and generate two files:  "blinding_key.csv" with the original and blinded file names, and "scoresheet.csv" with just the blinded filenames.  Add columns as needed to scoresheet, for example, runx_count.  Then run bb_unblind to rejoin scoresheet to the key and generate an unblinded result file.
#' @param analysis_file The analysis file for the experiment.  It should contain 1 line for every biological sample and should have a filename for every file to be blinded.
#' @param file_column The column name in the analysis_file with the files to be blinded.
#' @param output_dir The linux-style file path for the directory that will hold the blinded images.  The directory will be created by the function.
#' @return nothing
#' @export
#' @importFrom stringr str_replace_all str_extract str_replace str_sub
#' @importFrom fs path dir_create path_file file_copy
#' @importFrom dplyr mutate across select pull filter arrange
#' @importFrom tidyr contains
#' @importFrom digest digest
#' @importFrom tibble tibble
#' @importFrom readr write_csv
bb_blind_images <- function(analysis_file, file_column, output_dir) {
  ts <- stringr::str_replace_all(Sys.time(), "[:punct:]|[:alpha:]|[:space:]", "")
  output_dir <- fs::path(paste0(output_dir, "_", ts))
  fs::dir_create(path = output_dir)
  analysis_file <-
    analysis_file |>
    dplyr::mutate(dplyr::across(tidyr::contains(file_column), bb_fix_file_path))
  filepaths <-
    analysis_file |>
    dplyr::select(filepaths = tidyr::contains(file_column)) |>
    dplyr::mutate(filepaths = bb_fix_file_path(filepaths)) |>
    dplyr::pull(filepaths)

  blinding_key <- map_dfr(
    .x = filepaths,
    .f = function(x, out = output_dir) {
      filename <- fs::path_file(x)
      extension <- stringr::str_extract(x, ".tif|.tiff|.jpeg|.jpg|.png")
      filestem <- stringr::str_replace(filename, ".tif|.tiff|.jpeg|.jpg|.png", "")
      hash <- digest::digest(x, algo = "md5", file = TRUE)
      hash_short <- stringr::str_sub(hash, end = 6)
      hash_int <- strtoi(paste0("0x", hash_short))
      hash_mod <- hash_int %% nrow(wordhash)
      hash_final <- wordhash |>
        dplyr::filter(index == hash_mod) |>
        dplyr::pull(word) #%>%
        # paste(., "_", hash_short)
      key <-
        tibble::tibble(source_file = x,
               blinded_file = paste0(out, "/", hash_final, extension))
      return(key)
    }
  )
  fs::file_copy(path = blinding_key$source_file, new_path = blinding_key$blinded_file, overwrite = FALSE)
  scoresheet <-
    blinding_key |>
    dplyr::select(-source_file) |>
    dplyr::arrange(blinded_file)
  readr::write_csv(blinding_key, file = fs::path(output_dir, "blinding_key.csv"))
  readr::write_csv(scoresheet, file = fs::path(output_dir, "scoresheet.csv"))
}
