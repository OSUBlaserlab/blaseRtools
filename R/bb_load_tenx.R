#' Load 10X Data Into CDS
#'
#' @param targz_file A character string of the file path to the multi pipestance directory
#' @param umi_cutoff Don't import cells with fewer UMIs than this value.  Defaults to 100.
#' @param sample_metadata_tbl A tibble in wide format with one line.  Col names indicate metadata variables to add.
#' @return A cell data set object.
#' @export
#' @import monocle3 tidyverse SingleCellExperiment
bb_load_tenx_targz <- function (targz_file,
                                  umi_cutoff = 100,
                                  sample_metadata_tbl = NULL) {
  if (!file.exists(targz_file))
    stop(
      "Could not find the .tar.gz file: '",
      targz_file,
      "'.\n         Please double-check if the file exists.\n"
    )
  # create a temporary working directory in the active project
  dir.create("temp")
  # copy the tarball
  cmd <-
    paste0("cp ", targz_file, " temp")
  message(cmd, "\n")
  system(cmd)
  targz <- list.files("temp", full.names = T)
  #unzip
  cmd <-
    paste0("tar -xvf ", targz,  " -C temp")
  message(cmd, "\n")
  system(cmd)
  features.loc <- "temp/features.tsv.gz"
  barcode.loc <- "temp/barcodes.tsv.gz"
  matrix.loc <- "temp/matrix.mtx.gz"
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
  unlink("temp", recursive = TRUE)
  return(gbm)
}


#' @title Load a 10X Genomics H5 File and Return a CDS
#' @description This function reads a 10X Genomics H5 file and returns a cell_data_set or CDS.  Important notes:  This is tested and should work for all single-genome and citeseq data sets.  For multigenome data, as long as the features are all contained in the same matrix and identified by a composite reference/gene identifier, it should also work.  In this case, the CDS will have to be filtered post-hoc using the sample_barcodes.csv data to get the appropriate species of cell.  Also:  this function takes in a specific H5 file from a unique biological sample.  So it should be wrapped in another function to map across all the samples in a dataset.  The wrapper needs to find the appropriate H5 file, e.g. filtered_feature_bc_matrix.h5 for files processed with cellranger count or sample_filtered_feature_bc_matrix.h5 for files processed using cellranger multi. This may change based on the cellranger version used.
#' @param filename Path to the h5 file.
#' @return A cell data set.
#' @rdname bb_load_tenx_h5
#' @export
#' @importFrom cli cli_abort
#' @importFrom hdf5r H5File existsGroup
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat as.sparse
#' @importFrom monocle3 new_cell_data_set
#' @importFrom SingleCellExperiment mainExpName
bb_load_tenx_h5 <- function (filename) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    cli::cli_abort("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    cli::cli_abort("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, "matrix")) {
    feature_slot <- "features/id"
    gene_names <- "features/name"
  } else {
    feature_slot <- "genes"
    gene_names <- "gene_names"
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    features <- infile[[paste0(genome, "/", feature_slot)]][]
    gene_names <- infile[[paste0(genome, "/", gene_names)]][]
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- Matrix::sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      repr = "T",
    )
    features <- make.unique(names = features, sep = "__")
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- Seurat::as.sparse(x = sparse.mat)
    if (infile$exists(name = paste0(genome, "/features"))) {
      types <- infile[[paste0(genome, "/features/feature_type")]][]
    }
    barcodes <-
      data.frame(barcode = barcodes[], row.names = barcodes[])
    features <-
      data.frame(
        id = features,
        gene_short_name = gene_names,
        data_type = types,
        row.names = features
      )
    output[[genome]] <-
      list(mat = sparse.mat,
           features = features,
           barcodes = barcodes)
  }
  infile$close_all()

  cds <-
    monocle3::new_cell_data_set(
      expression_data = output$matrix$mat,
      cell_metadata = output$matrix$barcodes,
      gene_metadata =  output$matrix$features
    )
  SingleCellExperiment::mainExpName(cds) <- "Gene Expression"

  if ("Antibody Capture" %in% types) {
    cds <- bb_split_citeseq(cds)
  }

  return(cds)
}
