#' Make a Copy of Image Files and Rename With File Hashes in Blinded Folder
#'
#' @description Will copy and rename  the files and generate two files:  "blinding_key.csv" with the original and blinded file names, and "scoresheet.csv" with just the blinded filenames.  Add columns as needed to scoresheet, for example, runx_count.  Then run bb_unblind to rejoin scoresheet to the key and generate an unblinded result file.
#' @param analysis_file The analysis file for the experiment.  It should contain 1 line for every biological sample and should have a filename for every file to be blinded.  It should be read from a csv using read_csv and in the form of a tibble.  If necessary it should be filtered upstream of this function to contain only the samples which you want to have blinded together.  Generally you would filter out unuseable images.  Also you would usually only wish to blind images from treatment groups from the same clutch.  So if you had 2 treatment groups replicated in 3 clutches, you would filter appropriately to generate 3 sets of blinded images across treatment groups.
#' @param file_column The column name in the analysis_file with the files to be blinded.
#' @param output_dir The linux-style file path for the directory that will hold the blinded images.  The directory will be created by the function.
#' @return nothing
#' @export
#' @import tidyverse digest fs
bb_blind_images <- function(analysis_file, file_column, output_dir) {
  ts <- str_replace_all(Sys.time(), "[:punct:]|[:alpha:]|[:space:]", "")
  output_dir <- paste0(output_dir, "_", ts)
  dir.create(output_dir)
  stopifnot("You need a column named useable in the analysis file." = "useable" %in% colnames(analysis_file))
  analysis_file <-
    analysis_file %>%
    mutate(across(contains(file_column), bb_fix_file_path))
  filepaths <-
    analysis_file %>%
    filter(useable) %>%
    select(filepaths = contains(file_column)) %>%
    mutate(filepaths = bb_fix_file_path(filepaths)) %>%
    pull(filepaths)

  blinding_key <- map_dfr(
    .x = filepaths,
    .f = function(x, out = output_dir) {
      filename <- fs::path_file(x)
      extension <- str_extract(x, ".tif|.tiff|.jpeg|.jpg|.png")
      filestem <- str_replace(filename, ".tif|.tiff|.jpeg|.jpg|.png", "")
      hash <- digest::digest(x, algo = "md5", file = TRUE)
      hash_short <- str_sub(hash, end = 6)
      hash_int <- strtoi(paste0("0x", hash_short))
      hash_mod <- hash_int %% nrow(wordhash)
      hash_final <- wordhash %>%
        filter(index == hash_mod) %>%
        pull(word) #%>%
        # paste(., "_", hash_short)
      key <-
        tibble(source_file = x,
               blinded_file = paste0(out, "/", hash_final, extension))
      return(key)
    }
  )
  fs::file_copy(path = blinding_key$source_file, new_path = blinding_key$blinded_file, overwrite = FALSE)
  scoresheet <-
    blinding_key %>% select(-source_file) %>% arrange(blinded_file)
  write_csv(blinding_key, file = paste0(output_dir, "/blinding_key.csv"))
  write_csv(scoresheet, file = paste0(output_dir, "/scoresheet.csv"))
}
