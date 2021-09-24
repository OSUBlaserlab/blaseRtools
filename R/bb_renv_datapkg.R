#' Install Or Update A Local R Data Package
#'
#' @param path Path to a directory containing one or more versions of the same data package.  If a directory is specified, this function will compare the currently installed version to the latest availaible version in the directory.  If there is a newer version available (based on version number), it will install this version.  If a binary is specifically requested it will install that one.
#' @return Returns nothing. Using renv, it installs the latest version of the binary datapackage or updates it if already installed.
#' @export
#' @import tidyverse renv
bb_renv_datapkg <- function(path) {
  if (str_detect(string = path, pattern = ".tar.gz")) {
    message(str_glue("Installing {path}.  There may be newer versions available."))
    renv::install(path)
  } else {
    latest_version <- file.info(list.files(path, full.names = T)) %>%
      as_tibble(rownames = "file") %>%
      filter(str_detect(file, pattern = ".tar.gz")) %>%
      arrange(desc(mtime)) %>%
      dplyr::slice(1) %>%
      pull(file) %>%
      str_split(pattern = "/") %>%
      map(tail, n = 1) %>%
      unlist()
    datapackage_stem <- str_replace(latest_version, "_.*", "")
    latest_version_number <- str_replace(latest_version, "^.*_", "")
    latest_version_number <-
      str_replace(latest_version_number, ".tar.gz", "")
    if (packageVersion(datapackage_stem) < latest_version_number) {
      message(str_glue("A newer data package version is available.  Installing {latest_version}."))
      if (str_sub(path, -1) == "/") {
        renv::install(paste0(path, latest_version))
      } else {
        renv::install(paste0(path, "/", latest_version))
      }

    } else {
      message(str_glue("Your current version of {datapackage_stem} is up to date."))
    }
  }
}
