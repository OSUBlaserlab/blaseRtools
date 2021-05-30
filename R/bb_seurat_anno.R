#' A function to generate automated cell labelings with Seurat 
# 
#' @param cds A cell data set object
#' @param reference Seurat reference data. 
#' @return A modified cds with Seurat cell assignments. 
#' @export
#' @import tidyverse Seurat monocle3 broom Matrix SingleCellExperiment SeuratDisk
#' @examples
bb_seurat_anno <- function(cds, reference) {
  seurat_reference <- LoadH5Seurat(reference)
  
  exprs_renamed_long <- broom::tidy(exprs(cds)) %>%
    as_tibble() %>%
    left_join(., rowData(cds) %>%
                as_tibble(),
              by = c("row" = "id")) %>%
    select(gene_short_name, column, value) %>%
    select(row = gene_short_name, column, value)
  exprs_renamed_df <- as.data.frame(exprs_renamed_long)
  exprs_renamed_rows <-
    exprs_renamed_df %>%
    select(row) %>%
    unique() %>%
    data.frame()
  exprs_renamed_cols <-
    exprs_renamed_df %>%
    select(column) %>%
    unique() %>%
    data.frame()
  exprs_renamed_df$rowIdx <-
    match(exprs_renamed_df$row,
          exprs_renamed_rows$row)
  exprs_renamed_df$colIdx <-
    match(exprs_renamed_df$column,
          exprs_renamed_cols$column)
  exprs_spars <- sparseMatrix(
    i = exprs_renamed_df$rowIdx,
    j = exprs_renamed_df$colIdx,
    x = as.integer(exprs_renamed_df$value),
    dims = c(nrow(exprs_renamed_rows), nrow(exprs_renamed_cols)),
    dimnames = list(exprs_renamed_rows$row, exprs_renamed_cols$column)
  )
  seurat_cds <-
    CreateSeuratObject(counts = exprs_spars, assay = "RNA")
  
  seurat_cds <- SCTransform(seurat_cds)
  
  anchors <- FindTransferAnchors(
    reference = seurat_reference,
    query = seurat_cds,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50,
    features = intersect(
      rownames(x = seurat_reference),
      VariableFeatures(object = seurat_cds)
    ),
    reference.assay = "SCT",
    query.assay = "SCT",
    verbose = T,
    recompute.residuals = F
  )
  
  seurat_cds <- TransferData(
    anchorset = anchors,
    reference = seurat_reference,
    query = seurat_cds,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    )
  )
  
  seurat_cds <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = seurat_reference,
    query = seurat_cds,
    new.reduction.name = "ref.spca"
  )
  
  seurat_cds <- ProjectUMAP(
    query = seurat_cds,
    query.reduction = "ref.spca",
    reference = seurat_reference,
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  
  seurat_metadata <-
    as_tibble(seurat_cds@meta.data, rownames = "barcode")
  seurat_dims <-
    as_tibble(seurat_cds@reductions[["ref.umap"]][[1:dim(seurat_cds@reductions[["ref.umap"]])[1]]], rownames = "barcode")
  seurat_data <- left_join(seurat_dims, seurat_metadata)
  cds_data <- left_join(colData(cds) %>% as_tibble(rownames = "barcode_sample"),
          seurat_data,
          by = c("barcode_sample" = "barcode"))
  colData(cds)$seurat_dim1 <- cds_data$refUMAP_1
  colData(cds)$seurat_dim2 <- cds_data$refUMAP_2
  colData(cds)$seurat_celltype_l1 <- cds_data$predicted.celltype.l1
  colData(cds)$seurat_celltype_l2 <- cds_data$predicted.celltype.l2
  return(cds) 
  
}
