#' Install Or Update A Local R Data Package
#'
#' @param path Path to a directory containing one or more versions of the same data package.  There must be only 1 set of packages versions in the directory e.g. "my.package_1.0.tar.gz", "my.package_1.1.tar.gz", etc.
#' @return Returns nothing. Using renv, it installs the latest version of the binary datapackage or updates it if already installed.
#' @export
#' @import tidyverse renv
bb_renv_datapkg <- function(path) {
  latest_version <- file.info(list.files(path, full.names = T)) %>%
    as_tibble(rownames = "file") %>%
    filter(str_detect(file, pattern = ".tar.gz")) %>%
    arrange(desc(mtime)) %>%
    dplyr::slice(1) %>%
    pull(file) %>%
    str_split(pattern = "/") %>%
    map(tail, n=1) %>%
    unlist()
  if (str_sub(path, -1) == "/") {
    renv::install(paste0(path, latest_version))
  } else {
    renv::install(paste0(path, "/", latest_version))
  }
}




