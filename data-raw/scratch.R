blaseRtemplates::project_data(path = "~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/datapkg")
bb_cellmeta(cds_WT_AML_bALL) |> glimpse()
theme_set(theme_cowplot())
bb_plot_genes_in_pseudotime(
  cds = filter_cds(
    cds_WT_AML_bALL,
    cells = bb_cellmeta(cds_WT_AML_bALL) |> filter(leiden_best_consensus_all == "tusi_MPP")
  ),
  gene_or_genes = c("Kit", "Mpo"),
  pseudotime_dim = "pseudotime_leiden_best_consensus_all_tusi_MPP", legend_title = "pseudotime", color_cells_by = "louvain"
) + panel_border()


