#' Use doubletfinder to model and mark doublets
#'
#' @param cds A cell data set object
#' @param doublet_prediction Predicted proportion of doublets fom 0 to 1
#' @param qc_table A table of qc calls from the blaseRtools qc function
#' @return A tibble of low- and high-confidence doublet calls by barcode
#' @export
#' @import tidyverse Seurat
bb_doubletfinder <-
  function(cds,
           doublet_prediction,
           qc_table,
           ncores = 1
           ) {
    # system(paste0("gunzip -k ", directory, "/*"))
    keepers <- qc_table |> filter(qc.any == FALSE) |> pull(barcode)
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
      param_sweep(seu_object, PCs = 1:10, sct = FALSE, num.cores = ncores)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.stats)

    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- seu_object@meta.data$seurat_clusters
    homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
    nExp_poi <- round(doublet_prediction * length(annotations))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seu_object_lowconf <-
      doublet_finder(
        seu_object,
        PCs = 1:10,
        pN = 0.25,
        pK = 0.09,
        nExp = nExp_poi,
        reuse.pANN = FALSE,
        sct = FALSE

      )

    seu_object_highconf <-
      doublet_finder(
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


#' @import Seurat
#' @import fields
#' @importFrom SeuratObject Version
#' @importFrom DoubletFinder parallel_paramSweep_v3
param_sweep <- function (seu,
                         PCs = 1:10,
                         sct = FALSE,
                         num.cores = 1) {
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <-
    round(nrow(seu@meta.data) / (1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data),
                                                 10000, replace = FALSE)]
    if (SeuratObject::Version(seu) >= "5.0.0") {
      data <- seu@assays$RNA$counts[, real.cells]
    } else {
      data <- seu@assays$RNA@counts[, real.cells]
    }
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)

    if (SeuratObject::Version(seu) >= "5.0.0") {
      data <- seu@assays$RNA$counts[, real.cells]
    } else {
      data <- seu@assays$RNA@counts[, real.cells]
    }
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <-
      mclapply(
        as.list(1:length(pN)),
        FUN = DoubletFinder::parallel_paramSweep_v3,
        n.real.cells,
        real.cells,
        pK,
        pN,
        data,
        orig.commands,
        PCs,
        sct,
        mc.cores = num.cores
      )
    stopCluster(cl)
  }
  else {
    output2 <-
      lapply(
        as.list(1:length(pN)),
        FUN = DoubletFinder::parallel_paramSweep_v3,
        n.real.cells,
        real.cells,
        pK,
        pN,
        data,
        orig.commands,
        PCs,
        sct
      )
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK,
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

#' @importFrom SeuratObject Version
#' @import fields Seurat KernSmooth
doublet_finder <- function(seu,
          PCs,
          pN = 0.25,
          pK,
          nExp,
          reuse.pANN = FALSE,
          sct = FALSE,
          annotations = NULL) {
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <-
      "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp,
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    if (SeuratObject::Version(seu) >= "5.0.0") {
      data <- seu@assays$RNA$counts[, real.cells]
    } else {
      data <- seu@assays$RNA@counts[, real.cells]
    }
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells / (1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...",
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- base::cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <-
        NormalizeData(
          seu_wdoublets,
          normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
          scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
          margin = orig.commands$NormalizeData.RNA@params$margin
        )
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(
        seu_wdoublets,
        selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
        loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
        clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
        mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
        dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
        num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
        binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
        nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
        mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
        dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff
      )
      print("Scaling data...")
      seu_wdoublets <-
        ScaleData(
          seu_wdoublets,
          features = orig.commands$ScaleData.RNA$features,
          model.use = orig.commands$ScaleData.RNA$model.use,
          do.scale = orig.commands$ScaleData.RNA$do.scale,
          do.center = orig.commands$ScaleData.RNA$do.center,
          scale.max = orig.commands$ScaleData.RNA$scale.max,
          block.size = orig.commands$ScaleData.RNA$block.size,
          min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block
        )
      print("Running PCA...")
      seu_wdoublets <-
        RunPCA(
          seu_wdoublets,
          features = orig.commands$ScaleData.RNA$features,
          npcs = length(PCs),
          rev.pca = orig.commands$RunPCA.RNA$rev.pca,
          weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
          verbose = FALSE
        )
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells,
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells,
                                             ncol = length(levels(
                                               doublet_types1
                                             ))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells)) / k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbors_that_are_doublets = neighbors[neighbors >
                                                    n_real.cells]
          if (length(neighbors_that_are_doublets) > 0) {
            neighbor_types[i,] <-
              table(doublet_types1[neighbors_that_are_doublets -
                                     n_real.cells]) + table(doublet_types2[neighbors_that_are_doublets -
                                                                             n_real.cells])
            neighbor_types[i,] <- neighbor_types[i,] / sum(neighbor_types[i,])
          }
          else {
            neighbor_types[i,] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <-
      "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <-
      pANN[rownames(seu@meta.data),
           1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp,
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors",
                              pN, pK, nExp, ct, sep = "_")] <-
          neighbor_types[,
                         ct]
      }
    }
    return(seu)
  }
}
