#' Make a Copy of Image Files and Rename With File Hashes in Blinded Folder
#'
#' @description Will copy and rename  the files and generate two files:  "blinding_key.csv" with the original and blinded file names, and "scoresheet.csv" with just the blinded filenames.  Add columns as needed to scoresheet, for example, runx_count.  Then run bb_unblind to rejoin scoresheet to the key and generate an unblinded result file.
#' @param source_dir The linux-style file path of the image files to be blinded.
#' @param output_dir The linux-style file path for the blinded images.  This will be created by the function.
#' @return nothing
#' @export
#' @import tidyverse digest
bb_blind_images <- function(source_dir, output_dir) {
  filepaths <-
    list.files(source_dir, full.names = T, pattern = ".tif|.tiff|.jpeg|.jpg|.png")
  filenames <-
    list.files(source_dir, full.names = F, pattern = ".tif|.tiff|.jpeg|.jpg|.png")
  dir.create(output_dir)
  blinding_key <- map2_dfr(
    .x = filepaths,
    .y = filenames,
    .f = function(x, y, src = source_dir, out = output_dir) {
      extension <- str_extract(y, ".tif|.tiff|.jpeg|.jpg|.png")
      hash <- digest::digest(x, algo = "md5", file = TRUE)
      hash <- str_sub(hash, end = 6)
      key <-
        tibble(source_file = y,
               blinded_file = paste0(hash, extension))
      file.copy(
        from = x,
        to = output_dir,
        recursive = F,
        overwrite = T,
        copy.mode = T,
        copy.date = T
      )
      file.rename(
        from = paste0(output_dir, "/", y),
        to = paste0(output_dir, "/", hash, extension)
      )
      return(key)
    }
  )
  scoresheet <-
    blinding_key %>% select(-source_file) %>% arrange(blinded_file)
  write_csv(blinding_key, file = paste0(output_dir, "/blinding_key.csv"))
  write_csv(scoresheet, file = paste0(output_dir, "/scoresheet.csv"))
}
