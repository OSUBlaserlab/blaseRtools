#' Load 10X Data Into CDS
#'
#' @param targz_file A character string of the file path to the multi pipestance directory
#' @param umi_cutoff Don't import cells with fewer UMIs than this value.  Defaults to 100.
#' @param sample_metadata_tbl A tibble in wide format with one line.  Col names indicate metadata variables to add.
#' @return A cell data set object.
#' @export
#' @import monocle3 tidyverse
bb_load_tenx_targz <- function (targz_file,
                                  umi_cutoff = 100,
                                  allowed_data_types = c("Gene Expression", "Antibody Capture"),
                                  sample_metadata_tbl = NULL) {
  if (!file.exists(targz_file))
    stop(
      "Could not find the .tar.gz file: '",
      targz_file,
      "'.\n         Please double-check if the file exists.\n"
    )
  # create a temporary working directory in the active project
  temp <- tempdir(check = T)
  # copy the tarball
  cmd <-
    paste0("cp ", targz_file, " ", temp)
  message(cmd, "\n")
  system(cmd)
  targz <- list.files(temp, full.names = T)
  #unzip
  cmd <-
    paste0("tar -xvf ", targz,  " -C ", temp)
  message(cmd, "\n")
  system(cmd)
  features.loc <- file.path(temp, "features.tsv.gz")
  barcode.loc <- file.path(temp, "barcodes.tsv.gz")
  matrix.loc <- file.path(temp, "matrix.mtx.gz")
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing")
  }
  if (!file.exists(features.loc)) {
    stop("Gene name or features file missing")
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing")
  }
  data <- Matrix::readMM(matrix.loc)
  feature.names = utils::read.delim(features.loc,
                                    header = FALSE,
                                    stringsAsFactors = FALSE)
  feature.names$V1 = make.unique(feature.names$V1)
  if (dim(data)[1] != length(feature.names[, 1])) {
    stop(sprintf(
      paste(
        "Mismatch dimension between gene file: \n\t %s\n and",
        "matrix file: \n\t %s\n"
      ),
      features.loc,
      matrix.loc
    ))
  }
  colnames(feature.names) = c("id", "gene_short_name", "data_type")
  rownames(data) = feature.names[, "id"]
  rownames(feature.names) = feature.names[, "id"]
  barcodes <-
    utils::read.delim(barcode.loc,
                      stringsAsFactors = FALSE,
                      header = FALSE)
  if (dim(data)[2] != length(barcodes[, 1])) {
    stop(sprintf(
      paste(
        "Mismatch dimension between barcode file: \n\t %s\n",
        "and matrix file: \n\t %s\n"
      ),
      barcode.loc,
      matrix.loc
    ))
  }
  barcodes$V1 = make.unique(barcodes$V1)
  colnames(data) = barcodes[, 1]
  pd = data.frame(barcode = barcodes[, 1], row.names = barcodes[,
                                                                1])
  data <- data[, Matrix::colSums(data) > umi_cutoff]
  pd <- pd[colnames(data), , drop = FALSE]
  gbm <-
    monocle3::new_cell_data_set(data, cell_metadata = pd, gene_metadata = feature.names)
  if (!is.null(sample_metadata_tbl)) {
    sample_metadata_tbl <- pivot_longer(sample_metadata_tbl, cols = everything())
    for (i in 1:nrow(sample_metadata_tbl)) {
      colData(gbm)$new <- sample_metadata_tbl$value[i]
      names(colData(gbm))[names(colData(gbm)) == "new"] <-
        sample_metadata_tbl$name[i]
    }
  }
  unlink(temp, recursive = TRUE)
  return(gbm)
}
