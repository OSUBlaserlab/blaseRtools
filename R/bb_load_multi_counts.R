#' Load 10X Multi data into a cds
#'
#' @param pipestance_path A character string of the file path to the multi pipestance directory
#' @param specimen The specimen name.
#' @param genome A character string for multi-genome samples
#' @param umi_cutoff Don't import cells with fewer UMIs than this value
#' @param allowed_data_types Control type of data to import.
#' @param sample_metadata_nvps A named vector or list of name+value pairs where name is the cell metadata column name to add and value is a sample-level metadata value.  These are applied to all cells in the sample (e.g.:  patient_id, treatment, timepoint).
#' @return A cell data set object.
#' @export
#' @import monocle3
bb_load_multi_counts <- function (pipestance_path = NULL,
            specimen = NULL,
            genome = NULL,
            umi_cutoff = 100,
            allowed_data_types = c("Gene Expression","Antibody Capture"),
	          sample_metadata_nvps = NULL

  ){
    if (!dir.exists(pipestance_path))
      stop(
        "Could not find the pipestance path: '",
        pipestance_path,
        "'.\n         Please double-check if the directory exists.\n"
      )
    od = file.path(pipestance_path, "outs/per_sample_outs", specimen)
    if (!dir.exists(od))
      stop(
        "Could not find the pipestance output directory: '",
        file.path(pipestance_path, "outs"),
        "'. Please double-check if the directory exists.\n"
      )
    v3p = file.path(od, "count/sample_feature_bc_matrix")
    v2p = file.path(od, "filtered_gene_bc_matrices")
    v3d = dir.exists(v3p)
    matrix_dir = ifelse(v3d, v3p, v2p)
    if (!dir.exists(matrix_dir))
      stop("Could not find directory: ", matrix_dir)
    if (v3d) {
      features.loc <- file.path(matrix_dir, "features.tsv.gz")
      barcode.loc <- file.path(matrix_dir, "barcodes.tsv.gz")
      matrix.loc <- file.path(matrix_dir, "matrix.mtx.gz")
      summary.loc <- file.path(od, "metrics_summary_csv.csv")
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
    if (v3d) {
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
    }
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
      new_cell_data_set(data, cell_metadata = pd, gene_metadata = feature.names)
   if (!is.null(sample_metadata_nvps)) {
     for (i in 1:length(sample_metadata_nvps)) {
		  colData(gbm)$new <- unname(sample_metadata_nvps[i])
		  names(colData(gbm))[names(colData(gbm)) == "new"] <- names(columns_to_add[i])
								  }
   }
    return(gbm)
  }
