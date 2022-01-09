#' Run pseudobulk analysis using Deseq2
#'
#' @param cds_deseq The cell data set object subset to analyze
#' @param replicate_variable The colData column holding the biological replicate variable identifiers:  one unique identifier per biological replicate.
#' @param class_variable  The colData column holding the class variable identifiers:  one unique identifier per class; two classes only
#' @return A list of results from pseudobulk analysis
#' @import tidyverse pheatmap Matrix.utils monocle3
bb_pseudobulk <- function(cds_deseq,
			  replicate_variable,
			  class_variable) {
  message("This function has been deprecated.  You should consider using bb_pseudobulk_mf for multifactorial pseudobulk analysis.")
  groups <-
    colData(cds_deseq) %>%
    as_tibble() %>%
    select(replicate = !!sym(replicate_variable), class = !!sym(class_variable))
  # get the aggregate counts
  aggregate_counts <-
    aggregate.Matrix(t(as.matrix(exprs(cds_deseq))), groupings = groups, fun = "sum")
  counts_matrix <- as.matrix(t(as.matrix(aggregate_counts)))
  # make the coldata for deseq
  samples <- colnames(counts_matrix)
  classes <- str_extract(colnames(counts_matrix), paste(groups$class %>% as.character() %>% unique(),collapse = "|"))
  coldata <- data.frame(classes, row.names = samples)
  stopifnot(all(rownames(coldata) == colnames(counts_matrix)))
  # make the deseq object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = coldata,
                                design = ~ classes)

  # do the thing
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  result <-
    as.data.frame(res)  %>% rownames_to_column(var = "id") %>% as_tibble() %>%
    left_join(., as_tibble(rowData(cds_deseq)[, c("id", "gene_short_name")]))

  #qc
  rld <- DESeq2::rlog(dds,blind = T)
  pca_plot <- DESeq2::plotPCA(rld, intgroup = "classes")

  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)

  # Plot heatmap
  heatmap <- pheatmap(rld_cor, annotation = coldata[, c("classes"), drop = F])
  dispersion <- DESeq2::plotDispEsts(dds)

  return_list <- list(res@elementMetadata@listData[["description"]],
                      result,
                      pca_plot,
                      heatmap,
                      dispersion
                      )
  return(return_list)
}
