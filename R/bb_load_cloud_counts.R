#' Load 10X Cloud Pipestances Into CDS
#'
#' @param pipestance_path A character string of the file path to the multi pipestance directory
#' @param specimen The biological specimen name.
#' @param genome A character string for multi-genome samples; currently not working as originally intended so leave as null.
#' @param umi_cutoff Don't import cells with fewer UMIs than this value
#' @param allowed_data_types Control type of data to import.
#' @param sample_metadata_tbl A tibble in wide format with one line.  Col names indicate metadata variables to add.
#' @return A cell data set object.
#' @export
#' @import monocle3 tidyverse
bb_load_cloud_counts <- function (pipestance_path = NULL,
                                  specimen = NULL,
                                  genome = NULL,
                                  umi_cutoff = 100,
                                  allowed_data_types = c("Gene Expression", "Antibody Capture"),
                                  sample_metadata_tbl = NULL) {
  if (!dir.exists(pipestance_path))
    stop(
      "Could not find the pipestance path: '",
      pipestance_path,
      "'.\n         Please double-check if the directory exists.\n"
    )
  od = file.path(pipestance_path, "per_sample_outs", specimen)
  if (!dir.exists(od))
    stop(
      "Could not find the pipestance output directory: '",
      file.path(pipestance_path, "outs"),
      "'. Please double-check if the directory exists.\n"
    )
  # create a temporary working directory in the active project
  dir.create("temp")
  # copy the tarball from the 10X directory
  cmd <-
    paste0("cp ", od, "/count/sample_feature_bc_matrix.tar.gz temp")
  message(cmd, "\n")
  system(cmd)
  #unzip
  cmd <-
    paste0("tar -xvf temp/sample_feature_bc_matrix.tar.gz -C temp")
  message(cmd, "\n")
  system(cmd)
  if (is.null(genome)) {
    features.loc <- "temp/features.tsv.gz"
    barcode.loc <- "temp/barcodes.tsv.gz"
    matrix.loc <- "temp/matrix.mtx.gz"
    summary.loc <- file.path(od, "metrics_summary.csv")
  }
  else {
    genome = get_genome_in_matrix_path(matrix_dir, genome)
    barcode.loc <- file.path(matrix_dir, genome, "barcodes.tsv")
    features.loc <- file.path(matrix_dir, genome, "genes.tsv")
    matrix.loc <- file.path(matrix_dir, genome, "matrix.mtx")
    summary.loc <- file.path(od, "metrics_summary.csv")
  }
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
  data_types = factor(feature.names$V3)
  allowed = data_types %in% (allowed_data_types)
  if (!is.null(genome)) {
    gfilter = grepl(genome, feature.names$V1)
    if (any(gfilter)) {
      allowed = allowed & grepl(genome, feature.names$V1)
    }
    else {
      message(
        paste(
          "Data does not appear to be from a multi-genome sample,",
          "simply returning all gene feature data without",
          "filtering by genome."
        )
      )
    }
  }
  data = data[allowed,]
  feature.names = feature.names[allowed, 1:2]
  colnames(feature.names) = c("id", "gene_short_name")
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
  unlink("temp", recursive = T)
  return(gbm)
}
