#' Rejoin Blinded Scores to Original File Names
#'
#' @description Will rejoin scoresheet and blinded key to produce unblinded results.  If you change the names of either of those files, they have to be provided as arguments to the function.  Otherwise keyfile and scorefile are optional.
#' @param directory The linux-style filepath of the folder containing the scoresheet and blinded key.
#' @param keyfile Optional:  filename of the key file.  Defaults to "blinding_key.csv".
#' @param scorefile Optional:  filename of the score file.  Defaults to "scoresheet.csv".
#' @return nothing
#' @export
#' @import tidyverse
bb_unblind_images <-
  function(directory,
           keyfile = "blinding_key.csv",
           scorefile = "scoresheet.csv") {
    key <- read_csv(paste0(directory, "/", keyfile))
    scores <- read_csv(paste0(directory, "/", scorefile))
    res <- left_join(key,
                     scores) %>% arrange(source_file)
    write_csv(res, file = paste0(directory, "/unblinded_result.csv"))
  }
