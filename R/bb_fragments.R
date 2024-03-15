#' @title Replace The Path to a Fragment File In a Signac Object
#' @description Often you will wish to package a Signac object or share the object to someone who does not have access to the same directories as you.  This is a problem because one of the internal sub-objects, the Fragment object, holds a file path to a directory containing an atac fragments file and its .tbi intex.  Often this file will only be accessbile to you.  This function allows you to replace the fragments file path that is held internally within the Signac object.  Just copy the files to a shared location and provide that file location to the function.
#'
#' For Signac objects with multiple internal Fragments objects, provide vector of file paths.
#'
#' For packaged fragments files, use system.file or fs::path_package to access this, usually from extdata.
#'
#' @param obj A signac object.
#' @param new_paths A character vector of new file paths.  You must have the same number of new file paths as Fragments sub-objects within the signac object.
#' @return A modified Signac object.
#' @seealso
#'  \code{\link[Signac]{c("Cells.Fragment", "CountFragments", "CreateFragmentObject", "FilterCells", "Fragment-class", "Fragments", "Fragments", "SplitFragments", "UpdatePath", "ValidateCells", "ValidateFragments", "ValidateHash", "head.Fragment", "subset.Fragment")}}, \code{\link[Signac]{CreateFragmentObject}}
#'  \code{\link[cli]{cli_abort}}
#'  \code{\link[fs]{file_access}}
#'  \code{\link[purrr]{map2}}
#' @rdname bb_fragment_replacement
#' @export
#' @importFrom Signac Fragments CreateFragmentObject
#' @importFrom cli cli_abort
#' @importFrom fs file_access
#' @importFrom purrr map2
bb_fragment_replacement <- function(obj,
                         new_paths) {
  n_new_paths <- length(new_paths)
  old <- Signac::Fragments(obj)
  n_old <- length(old)
  if (n_new_paths != n_old) {
    cli::cli_abort("The number of new paths must match the number of fragment objects in obj.")
  }

  if (!all(fs::file_access(new_paths))) {
    cli::cli_abort("You do not have access to this fragments file.")
  }

  if (n_new_paths == 1) {
    cellnames <- old@cells
    new <- Signac::CreateFragmentObject(path = new_paths,
                                        cells = cellnames,
                                        validate.fragments = TRUE)
  } else {
    # replace more than one fragments file as a list
    new <- purrr::map2(.x = old,
                .y = new_paths,
                .f = \(x, y) {
                  cellnames <- x@cells
                  new_frag <- Signac::CreateFragmentObject(path = y,
                                                           cells = cellnames,
                                                           validate.fragments = TRUE)
                })

  }

  Signac::Fragments(obj) <- NULL
  Signac::Fragments(obj) <- new

  return(obj)
}
