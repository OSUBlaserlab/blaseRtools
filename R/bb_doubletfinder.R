#' Use doubletfinder to modela and mark doublets
#' 
#' @param cds A cell data set object
#' @param doublet_prediction Predicted proportion of doublets fom 0 to 1
#' @param qc_table A table of qc calls from the blaseRtools qc function
#' @return A tibble of low- and high-confidence doublet calls by barcode 
#' @export
#' @import tidyverse DoubletFinder Seurat
#' @examples
bb_doubletfinder <-
  function(cds,
           doublet_prediction,
           qc_table
           ) {
    # system(paste0("gunzip -k ", directory, "/*"))
    keepers <- qc_table %>% filter(qc.any == FALSE) %>% pull(barcode)
    #seu_data <- Read10X(data.dir = directory)
    seu_object <- CreateSeuratObject(exprs(cds))
    seu_object@meta.data$barcode <- rownames(seu_object@meta.data)
    seu_object <- subset(seu_object, subset = barcode %in% keepers)
    seu_object <- NormalizeData(seu_object)
    seu_object <-
      FindVariableFeatures(seu_object,
                           selection.method = "vst",
                           nfeatures = 2000)
    all.genes  <- rownames(seu_object)
    seu_object <- ScaleData(seu_object, features = all.genes)
    seu_object <-
      RunPCA(seu_object, features = VariableFeatures(object = seu_object))
    seu_object <- FindNeighbors(seu_object)
    seu_object <- FindClusters(seu_object)
    seu_object <- RunUMAP(seu_object, dims = 1:10)
    
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list <-
      paramSweep_v3(seu_object, PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- seu_object@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(doublet_prediction * length(annotations))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seu_object_lowconf <-
      doubletFinder_v3(
        seu_object,
        PCs = 1:10,
        pN = 0.25,
        pK = 0.09,
        nExp = nExp_poi,
        reuse.pANN = FALSE,
        sct = FALSE
      )
    seu_object_highconf <-
      doubletFinder_v3(
        seu_object,
        PCs = 1:10,
        pN = 0.25,
        pK = 0.09,
        nExp = nExp_poi.adj,
        reuse.pANN = FALSE,
        sct = FALSE
      )
    
    ## extract the barcode and singlet/doublet calls
    solm <-
      as_tibble(seu_object_lowconf@meta.data) %>% select(barcode, doubletfinder_low_conf = contains("classifications"))
    sohm <-
      as_tibble(seu_object_highconf@meta.data) %>% select(barcode, doubletfinder_high_conf = contains("classifications"))
    
    ## join together and export as one
    
    return(full_join(solm, sohm))
  }
