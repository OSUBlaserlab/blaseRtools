devtools::load_all()
blaseRtemplates::project_data("~/network/X/Labs/Blaser/share/resources/datapkg/blaseRextras/")

test <- bb_monocle_anno(cds_qry = vignette_cds, cds_ref = pbmc_ref, labels = c("celltype.l1", "celltype.l2", "celltype.l3"))

