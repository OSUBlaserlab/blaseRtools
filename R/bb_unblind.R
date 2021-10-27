#' Rejoin Blinded Scores to Original File Names
#'
#' @description Will rejoin scoresheet and blinded key to produce unblinded results.  If you change the names of either of those files, they have to be provided as arguments to the function.  Otherwise keyfile and scorefile are optional.
#' @param directory The linux-style filepath of the folder containing the scoresheet and blinded key.
#' @param keyfile Optional:  filename of the key file.  Defaults to "blinding_key.csv".
#' @param scorefile Optional:  filename of the score file.  Defaults to "scoresheet.csv".
#' @param analysis_file Optional:  complete file path to the the unblinded analysis sheet.  If NULL (default), the function will only produce unblinded_result.csv which you will have to manually copy to the main analysis file.  If a valid file is specified, it will left_join unblinded results to the analysis_file.  In the process, it will necessarily convert windows file paths to linux-style file paths.  Samples not included in the blinding should return with NA values for the added columns.  Since this will overwrite a data file always commit your analysis_file before running with this specification.  The function will ask you to commit before writing.
#' @param file_column Optional:  The column in analysis_file containing file paths for the files that were blinded.
#' @return nothing
#' @export
#' @import tidyverse
bb_unblind_images <-
  function(directory,
           keyfile = "blinding_key.csv",
           scorefile = "scoresheet.csv",
           analysis_file = NULL,
           file_column = NULL) {
    key <- read_csv(paste0(directory, "/", keyfile))
    scores <- read_csv(paste0(directory, "/", scorefile))
    res <- left_join(key,
                     scores) %>%
      arrange(source_file)
    write_csv(res, file = paste0(directory, "/unblinded_result.csv"))
    if (!is.null(analysis_file)) {
      response <-
        readline(
          "This function will write to an existing data file.  Have you previously committed the analysis_file (y/n)? Unless y is returned the function will stop.   \n"
        )
      stopifnot(
        "Commit the analysis_file before running this function with analysis_file specified." = response == "y"
      )
      stopifnot("The analysis file is not present at this location." = file.exists(analysis_file))
      stopifnot("You must specify a column with the source file location." = !is.null(file_column))
      analysis_file <-
        read_csv(analysis_file)  %>%
        mutate(across(contains(file_column), bb_fix_file_path)) %>%
        rename(source_file = contains(file_column))
      analysis_joined <-
        left_join(analysis_file, res)
      write_csv(analysis_joined, file = analysis_file)

    }
  }
