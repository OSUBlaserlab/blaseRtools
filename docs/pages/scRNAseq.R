## ---- include = FALSE------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)


## ----setup, results = "hide"-----------------------------------------------------------------
# Attach the packages you will need for the analysis.
library(blaseRtools)
library(blaseRdata)
library(monocle3)
library(Seurat)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size = 12))


## --------------------------------------------------------------------------------------------
# Read in and inspect the configuration file.
vignette_config <- read_csv(system.file("extdata/vignette_config.csv", 
                                        package = "blaseRdata"), 
                            col_type = list(.default = col_character()))
vignette_config



## --------------------------------------------------------------------------------------------
# Fix the windows-style file path.
vignette_config <- vignette_config %>%
  mutate(targz_path = bb_fix_file_path(targz_path))
vignette_config


## ---- eval = FALSE---------------------------------------------------------------------------
## # Generate a list of CDS objects using purrr::map
## cds_list <- map(
##   .x = vignette_config$sample,
##   .f = function(x, conf = vignette_config) {
##     conf_filtered <- conf %>%
##       filter(sample == x)
##     cds <- bb_load_tenx_targz(targz_file = conf_filtered$targz_path,
##                               sample_metadata_tbl = conf_filtered %>%
##                                 select(-c(sample, targz_path))
##                               )
##     return(cds)
##   }
## ) %>%
##   set_names(nm = vignette_config$sample)


## ----eval = FALSE----------------------------------------------------------------------------
## # generate a list of qc results for individual CDS objects
## vig_qc_res <- pmap(.l = list(cds = cds_list,
##                              cds_name = names(cds_list),
##                              genome = rep("human", times = length(cds_list))),
##                    .f = bb_qc
##                    ) %>%
##   set_names(nm = names(cds_list))
## 


## ---- eval = FALSE, dev='png', fig.align='center', fig.width=4.5, fig.height=3.5-------------
## vig_qc_res$chromium_controller[3]
## vig_qc_res$chromium_controller[4]


## ----eval = FALSE----------------------------------------------------------------------------
## 
## # gets the number of cells in each cds and divides it by 100000
## anticipated_doublet_rate <- unlist(map(cds_list, ncol))/100000
## 
## # extracts the first element of the qc result list for each cds
## qc_calls <- map(vig_qc_res,1)
## 
## # generates a list of tables with doubletfinder results
## doubletfinder_list <-
##   pmap(
##     .l = list(
##       cds = cds_list,
##       doublet_prediction = anticipated_doublet_rate,
##       qc_table = qc_calls
##     ),
##     .f = bb_doubletfinder
##   ) %>%
##   set_names(names(cds_list))
## 
## 


## ---- eval = FALSE---------------------------------------------------------------------------
## # rejoins doubletfinder and qc data onto the list of CDS objects
## cds_list <- pmap(
##   .l = list(
##     cds = cds_list,
##     qc_data = qc_calls,
##     doubletfinder_data = doubletfinder_list
##   ),
##   .f = bb_rejoin
## )


## ----eval = FALSE----------------------------------------------------------------------------
## # Merge the CDS list into a single CDS
## vignette_cds <- monocle3::combine_cds(cds_list = cds_list)


## ----eval = FALSE----------------------------------------------------------------------------
## # Remove mitochondrial and ribosomal genes.
## vignette_cds <-
##   vignette_cds[rowData(vignette_cds)$gene_short_name %notin% hg38_remove_genes,]


## ----eval = FALSE----------------------------------------------------------------------------
## # Remove the low-quality cells
## vignette_cds <- vignette_cds[,colData(vignette_cds)$qc.any == FALSE]
## 
## # Remove the high-confidence doublets
## vignette_cds <-
##   vignette_cds[,colData(vignette_cds)$doubletfinder_high_conf == "Singlet"]
## 
## # Now remove the qc and doubletfinder columns from the cell metadata because we are done with them.
## colData(vignette_cds)$qc.any <- NULL
## colData(vignette_cds)$doubletfinder_low_conf <- NULL
## colData(vignette_cds)$doubletfinder_high_conf <- NULL
## 


## ----eval = FALSE----------------------------------------------------------------------------
## # Calculate the PCA dimensions
## vignette_cds <- preprocess_cds(vignette_cds)


## ----eval = FALSE----------------------------------------------------------------------------
## # Calculate UMAP dimensions
## vignette_cds <- reduce_dimension(vignette_cds, cores = 40)


## --------------------------------------------------------------------------------------------
# Cell metadata
bb_cellmeta(vignette_cds) 

# Gene metadata
bb_rowmeta(vignette_cds)


## ----eval = FALSE----------------------------------------------------------------------------
## # Align samples according to the equipment metadata column
## vignette_cds_aligned_temp <- bb_align(vignette_cds, align_by = "sample")
## 
## # Replace the original CDS with the Aligned CDS
## vignette_cds <- vignette_cds_aligned_temp
## rm(vignette_cds_aligned_temp)


## ---- dev='png', fig.align='center', fig.width=6.0, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, var = "sample")


## --------------------------------------------------------------------------------------------
bb_cellmeta(vignette_cds)


## ---- dev='png', fig.align='center', fig.width=6.0, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, var = "sample", 
            alt_dim_x = "prealignment_dim1", 
            alt_dim_y = "prealignment_dim2")


## ----eval = FALSE----------------------------------------------------------------------------
## # Identify clusters and calculate top markers
## marker_file <- tempfile()
## vignette_cds <- bb_triplecluster(vignette_cds, n_top_markers = 50, outfile = marker_file, n_cores = 8)
## vignette_top_markers <- read_csv(marker_file)


## --------------------------------------------------------------------------------------------
vignette_top_markers


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, var = "partition")
bb_var_umap(vignette_cds, var = "leiden")
bb_var_umap(vignette_cds, var = "louvain")


## ----eval = FALSE----------------------------------------------------------------------------
## # Identify gene modules and add them to the gene metadata.
## vignette_cds <- bb_gene_modules(vignette_cds, n_cores = 24)


## --------------------------------------------------------------------------------------------
bb_rowmeta(vignette_cds)


## ----eval = FALSE----------------------------------------------------------------------------
## # Annotate the PBMC data
## vignette_cds <- bb_seurat_anno(vignette_cds, reference = "~/workspace_pipelines/sc_refdata/seurat_pbmc_reference_20210506/pbmc_multimodal.h5seurat")


## --------------------------------------------------------------------------------------------
bb_cellmeta(vignette_cds)


## ---- dev='png', fig.align='center', fig.width=5.5, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, 
            var = "seurat_celltype_l1", 
            alt_dim_x = "seurat_dim1", 
            alt_dim_y = "seurat_dim2", 
            overwrite_labels = TRUE, 
            group_label_size = 4)


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, var = "seurat_celltype_l1")


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, var = "leiden", plot_title = "Leiden Clusters")


## --------------------------------------------------------------------------------------------
leiden_seurat <- bb_cellmeta(vignette_cds) %>%
  group_by(leiden, seurat_celltype_l1) %>%
  summarise(n = n())
leiden_seurat


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
ggplot(leiden_seurat, 
       mapping = aes(x = leiden, 
                     y = n, 
                     fill = seurat_celltype_l1)) +
  geom_bar(position="fill", stat="identity") + 
  scale_fill_brewer(palette = "Set1")


## ----eval = FALSE----------------------------------------------------------------------------
## # Recode the leiden clusters with our cell assignments
## colData(vignette_cds)$leiden_assignment <- recode(colData(vignette_cds)$leiden,
##                                                   "1" = "T/NK",
##                                                   "2" = "DC/Mono",
##                                                   "3" = "B")


## ---- dev='png', fig.align='center', fig.width=6.0, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, var = "leiden_assignment")


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, 
            var = "leiden_assignment", 
            value_to_highlight = "T/NK", 
            cell_size = 2, 
            foreground_alpha = 0.4)


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, 
            var = "leiden_assignment", 
            value_to_highlight = "T/NK", 
            palette = "green4", 
            cell_size = 2, 
            foreground_alpha = 0.4)


## ---- dev='png', fig.align='center', fig.width=7.0, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, 
            var = "density", 
            facet_by = "equipment", 
            cell_size = 2, 
            plot_title = "Local Cell Density")


## ---- dev='png', fig.align='center', fig.width=7.0, fig.height=3.5---------------------------
bb_var_umap(vignette_cds, 
            var = "local_n", 
            nbin = 10, sample_equally = T, 
            facet_by = "equipment", 
            cell_size = 2, 
            plot_title = "Local N in 10 Bins")

bb_var_umap(vignette_cds, 
            var = "log_local_n", 
            nbin = 10, 
            sample_equally = T, 
            facet_by = "equipment", 
            cell_size = 2, 
            plot_title = "Log Local N in 10 Bins")


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_cluster_representation(cds = vignette_cds, 
                          cluster_var = "leiden_assignment", 
                          class_var = "equipment", 
                          experimental_class = "X", 
                          control_class = "chromium", 
                          return_value = "plot")


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_cluster_representation(cds = vignette_cds, 
                          cluster_var = "leiden_assignment", 
                          class_var = "equipment", 
                          experimental_class = "X", 
                          control_class = "chromium", 
                          return_value = "table")



## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_gene_umap(vignette_cds, 
             gene_or_genes = "CD3D")


## ---- dev='png', fig.align='center', fig.width=7.5, fig.height=3-----------------------------
bb_gene_umap(vignette_cds, gene_or_genes = c("CD19", "CD3D", "CD14"))



## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
agg_genes <- 
  bb_rowmeta(vignette_cds) %>%
  select(id, module_labeled) %>%
  filter(module_labeled == "Module 1")
  

bb_gene_umap(vignette_cds, 
             gene_or_genes = agg_genes)


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_gene_dotplot(vignette_cds, 
                markers = c("CD3E", "CD14", "CD19"), 
                group_cells_by = "leiden_assignment")


## ---- dev='png', fig.align='center', fig.width=4.5, fig.height=3.5---------------------------
bb_gene_dotplot(vignette_cds, 
                markers = c("CD3E", "CD14", "CD19"), 
                group_cells_by = "leiden_assignment", 
                gene_ordering = c("CD19", "CD14", "CD3E"),
                group_ordering = c("T/NK", "DC/Mono", "B"))



## ---- dev='png', fig.align='center', fig.width=6.0, fig.height=3.5---------------------------
cell_groupings <-
  tribble(
  ~aesthetic,~variable, ~value, ~level,
  "facet","equipment", "chromium", 1,
  "facet","equipment", "X", 2,
  "axis", "leiden_assignment", "B", 1,
  "axis", "leiden_assignment", "T/NK", 2,
  "axis", "leiden_assignment", "DC/Mono", 3 
  )
    

bb_gene_dotplot(vignette_cds, 
                markers = c("CD3E", "CD14", "CD19"), 
                group_cells_by = "multifactorial", 
                gene_ordering = c("CD19", "CD14", "CD3E"),
                group_ordering = cell_groupings)


## --------------------------------------------------------------------------------------------
vignette_exp_design <- 
  bb_cellmeta(vignette_cds) %>%
  group_by(sample, leiden_assignment) %>%
  summarise()
vignette_exp_design


## ---- eval=FALSE-----------------------------------------------------------------------------
## vignette_pseudobulk_res <-
##   bb_pseudobulk_mf(cds = vignette_cds,
##                    pseudosample_table = vignette_exp_design,
##                    design_formula = "~ leiden_assignment",
##                    result_recipe = c("leiden_assignment", "T/NK", "DC/Mono"))


## --------------------------------------------------------------------------------------------
# Detailed column headers for the results tables.
vignette_pseudobulk_res$Header 


## --------------------------------------------------------------------------------------------
# Differential expression results.  Positive L2FC indicates up in T/NK vs DC/Mono

vignette_pseudobulk_res$Result %>%
  filter(log2FoldChange > 0) %>%
  arrange(padj)
  
vignette_pseudobulk_res$Result %>%
  filter(log2FoldChange < 0) %>%
  arrange(padj)


## --------------------------------------------------------------------------------------------
vignette_regression_res <- bb_monocle_regression(cds = vignette_cds, 
                      gene_or_genes = c("CD19", "CD3D", "CD14"), 
                      form = "~leiden_assignment")
vignette_regression_res



## --------------------------------------------------------------------------------------------
vignette_regression_res$term

