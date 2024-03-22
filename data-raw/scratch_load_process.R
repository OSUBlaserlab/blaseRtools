# this is a general template for loading and processing 10X scrnaseq data
# modify as necessary
library("blaseRtools")
library("monocle3")
library("Seurat")
library("tidyverse")

# read in the config file -------------------------------------------------
analysis_configs <- read_csv(fs::path("~/network/X\\Labs\\Blaser\\share\\collaborators\\freud_nk\\pipestance\\caprice\\analysis_configs.csv"))
analysis_configs

# Generate a list of CDS objects using purrr::map ------------------------
cds_list <- map(.x = analysis_configs$sample[3],
                .f = \(x, conf = analysis_configs) {
                  conf_filtered <- conf |>
                    filter(sample == x)
                  h5 <- list.files(
                    conf_filtered$pipestance_path,
                    pattern = "filtered_feature_bc_matrix.h5",
                    recursive = T,
                    full.names = T
                  )
                  cds <- bb_load_tenx_h5(filename = h5,
                                         sample_metadata_tbl = conf_filtered |>
                                           select(-c(sample, pipestance_path)))
                  return(cds)
                })

test <- SingleCellExperiment::swapAltExp(cds_list[[1]], "Multiplexing Capture")
test@assays@data$counts
cds_list[[1]]@assays$ADT$counts[1:5,1:5]
cds_list[[1]]@assays$ADT@layers$data[1:5,1:5]
cds_list[[1]]@assays$ADT@layers$counts[1:5,1:5]

# generate a list of qc results for individual CDS objects ---------------
ind_qc_res <- pmap(.l = list(
  cds = cds_list,
  cds_name = names(cds_list),
  genome = rep("human", times = length(cds_list))
),
.f = bb_qc) %>%
  set_names(nm = names(cds_list))


# gets the number of cells in each cds and divides it by 100000 ---------
anticipated_doublet_rate <- unlist(map(cds_list, ncol)) / 100000

# extracts the first element of the qc result list for each cds ---------
qc_calls <- map(ind_qc_res, 1)

# generates a list of tables with doubletfinder results -----------------
doubletfinder_list <-
  pmap(
    .l = list(
      cds = cds_list,
      doublet_prediction = anticipated_doublet_rate,
      qc_table = qc_calls
    ),
    .f = bb_doubletfinder
  ) %>%
  set_names(names(cds_list))


# rejoins doubletfinder and qc data onto the list of CDS objects -------
cds_list_rejoined <- pmap(
  .l = list(
    cds = cds_list,
    qc_data = qc_calls,
    doubletfinder_data = doubletfinder_list
  ),
  .f = bb_rejoin
)

# Merge the CDS list into a single CDS ---------------------------------
cds_main <- monocle3::combine_cds(cds_list = cds_list_rejoined)

# Remove mitochondrial and ribosomal genes. ----------------------------
cds_main <-
  cds_main[rowData(cds_main)$gene_short_name %notin%
             blaseRdata::hg38_remove_genes, ]

# Remove the low-quality cells ----------------------------------------
cds_main <- cds_main[, colData(cds_main)$qc.any == FALSE]


# Uncomment and run the following section if you want to preview the CDS
# before removing predicted doublets and aligning.

# # option to preview cds ---------------------------------------------------
# # Calculate the PCA dimensions
# cds_main <- preprocess_cds(cds_main,
#                            use_genes = bb_rowmeta(cds_main) |>
#                              filter(data_type == "Gene Expression") |>
#                              pull(id))
# # Calculate UMAP dimensions
# cds_main <- reduce_dimension(cds_main, cores = 40)


# Remove the high-confidence doublets ------------------------------------
cds_main <-
  cds_main[, colData(cds_main)$doubletfinder_high_conf == "Singlet"]

# Remove the qc and doubletfinder columns from the cell metadata --------
colData(cds_main)$qc.any <- NULL
colData(cds_main)$doubletfinder_low_conf <- NULL
colData(cds_main)$doubletfinder_high_conf <- NULL

# Calculate the PCA dimensions ------------------------------------------
cds_main <- preprocess_cds(cds_main,
                           use_genes = bb_rowmeta(cds_main) |>
                             filter(data_type == "Gene Expression") |>
                             pull(id))
# Calculate UMAP dimensions ---------------------------------------------
cds_main <- reduce_dimension(cds_main, cores = 40)

# align by batch --------------------------------------------------------
# edit the align_by parameter accordingly
cds_main <- bb_align(cds_main, align_by = "batch")

# Identify clusters and calculate top markers ---------------------------
dir.create("data")
cds_main <-
  bb_triplecluster(
    cds_main,
    n_top_markers = 50,
    outfile = "data/cds_main_top_markers.csv",
    n_cores = 8
  )
cds_main_top_markers <- read_csv("data/cds_main_top_markers.csv")

# Identify gene modules and add them to the gene metadata. ---------------
cds_main <- bb_gene_modules(cds_main, n_cores = 24)

# align to seurat reference ---------------------------------------------
# NB:  human PBMC only
cds_main <-
  bb_seurat_anno(
    cds_main,
    reference = system.file("extdata", "pbmc_multimodal.h5seurat",
                            package = "blaseRextras")
  )

# save the objects -------------------------------------------------------
save(analysis_configs, file = "data/analysis_configs.rda", compress = "bzip2")
save(cds_main, file = "data/cds_main.rda", compress = "bzip2")
save(cds_main_top_markers, file = "data/cds_main_top_markers.rda", compress = "bzip2")
save(seq_qc, file = "data/seq_qc.rda", compress = "bzip2")
save(ind_qc_res, file = "data/ind_qc_res.rda", compress = "bzip2")
