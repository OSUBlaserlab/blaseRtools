#' A function to generate automated cell labelings with Seurat
#
#' @param cds A cell data set object
#' @param reference Seurat reference data.
#' @return A modified cds with Seurat cell assignments.
#' @export
#' @importFrom SeuratDisk LoadH5Seurat
#' @importFrom monocle3 exprs
#' @importFrom dplyr left_join pull
#' @importFrom tibble tibble as_tibble
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom Seurat SCTransform FindTransferAnchors VariableFeatures TransferData IntegrateEmbeddings ProjectUMAP
#' @importFrom SummarizedExperiment colData
bb_seurat_anno <- function(cds, reference) {
  seurat_reference <- SeuratDisk::LoadH5Seurat(reference)
  exprs <- monocle3::exprs(cds)
  rownames(exprs) <-
    dplyr::left_join(tibble::tibble(feature_id = rownames(exprs)),
                     bb_rowmeta(cds)) |>
    dplyr::pull(gene_short_name)
  seurat_cds <- SeuratObject::CreateSeuratObject(counts = exprs, assay = "RNA")
  seurat_cds <- Seurat::SCTransform(seurat_cds)
  anchors <- Seurat::FindTransferAnchors(
    reference = seurat_reference,
    query = seurat_cds,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50,
    features = intersect(
      rownames(x = seurat_reference),
      Seurat::VariableFeatures(object = seurat_cds)
    ),
    reference.assay = "SCT",
    query.assay = "SCT",
    verbose = T,
    recompute.residuals = F
  )
  seurat_cds <-
    Seurat::TransferData(
      anchorset = anchors,
      reference = seurat_reference,
      query = seurat_cds,
      refdata = list(
        celltype.l1 = "celltype.l1",
        celltype.l2 = "celltype.l2",
        predicted_ADT = "ADT"
      )
    )
  seurat_cds <-
    Seurat::IntegrateEmbeddings(
      anchorset = anchors,
      reference = seurat_reference,
      query = seurat_cds,
      new.reduction.name = "ref.spca"
    )
  seurat_cds <-
    Seurat::ProjectUMAP(
      query = seurat_cds,
      query.reduction = "ref.spca",
      reference = seurat_reference,
      reference.reduction = "spca",
      reduction.model = "wnn.umap"
    )
  seurat_metadata <-
    tibble::as_tibble(seurat_cds@meta.data, rownames = "barcode")
  seurat_dims <-
    tibble::as_tibble(seurat_cds@reductions[["ref.umap"]][[1:dim(seurat_cds@reductions[["ref.umap"]])[1]]],
              rownames = "barcode")
  seurat_data <- dplyr::left_join(seurat_dims, seurat_metadata)
  cds_data <-
    dplyr::left_join(
      SummarizedExperiment::colData(cds) |>
        tibble::as_tibble(rownames = "barcode_sample"),
      seurat_data,
      by = c(barcode_sample = "barcode")
    )
  SummarizedExperiment::colData(cds)$seurat_dim1 <- cds_data$refUMAP_1
  SummarizedExperiment::colData(cds)$seurat_dim2 <- cds_data$refUMAP_2
  SummarizedExperiment::colData(cds)$seurat_celltype_l1 <- cds_data$predicted.celltype.l1
  SummarizedExperiment::colData(cds)$seurat_celltype_l2 <- cds_data$predicted.celltype.l2
  return(cds)
}

