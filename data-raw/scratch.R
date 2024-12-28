devtools::load_all()
blaseRtemplates::project_data("~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/datapkg/")

cds_test <- filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> dplyr::filter(specimen == "cll1_baseline"))

options(future.globals.maxSize= 5242880000)
test_res <- bb_cellchat(cds_test, group_var = "seurat_l2_leiden_consensus", species = "human", n_cores = 4)

test <- bb_monocle_anno(cds_qry = vignette_cds, cds_ref = pbmc_ref, labels = c("celltype.l1", "celltype.l2", "celltype.l3"))

