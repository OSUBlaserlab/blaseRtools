bb_cellmeta(cds_main_stables_refiltered) |> glimpse()
bb_var_umap(cds_main_stables_refiltered, "density")
bb_var_umap(cds_main_stables_refiltered,
            "da_score_condition_control_hgfa",
            hexify = TRUE,
            n_hexbins = 10)
bb_var_umap(cds_main_stables_refiltered,
            "sample",
            hexify = TRUE)

p + scale_fill_gradient2()
View(p)
devtools::load_all()
test <- bb_cellchat_heatmap(cellchat_cds_main_stables_refiltered)

plot_grid(test)
test |>
  count(source_target, interaction_name_2) |> View()
test |> filter(interaction_name_2 == "si:ch73-22o12.1 - si:ch73-22o12.1") |> filter(source_target == "HSC/Mk -> HSC/Mk") |> View()
