

pseudotime_cds_1 <-
  bb_pseudotime(
    vignette_cds,
    cluster_variable = "leiden",
    cluster_value = "1",
    calculate_pseudotime = TRUE
  )
pseudotime_cds_1 <- bb_pseudotime(pseudotime_cds_1,
                                  cluster_variable = "leiden",
                                  cluster_value = "2")

bb_cellmeta(pseudotime_cds_1) |> glimpse()
bb_var_umap(pseudotime_cds_1, "leiden", overwrite_labels = TRUE)
bb_var_umap(pseudotime_cds_1,
            "pseudotime_leiden_2",
            show_trajectory_graph = TRUE,
            trajectory_graph_color = "red",
            trajectory_graph_segment_size = 2,
            label_principal_points = FALSE,
            label_roots = TRUE,
            label_leaves = FALSE,
            label_branch_points = TRUE)
bb_var_umap(filtered,
            "pseudotime_leiden_1",
            show_trajectory_graph = TRUE,
            trajectory_graph_color = "red",
            trajectory_graph_segment_size = 2,
            label_principal_points = FALSE,
            label_roots = TRUE,
            label_leaves = TRUE,
            label_branch_points = TRUE)
monocle3::plot_cells(pseudotime_cds_1, color_cells_by = "pseudotime", label_leaves = FALSE)


filtered <- filter_cds(pseudotime_cds_1, cells = bb_cellmeta(pseudotime_cds_1) |> filter(leiden == "1"))

bb_var_umap(filtered, "pseudotime_leiden_1", show_trajectory_graph = TRUE)

str(filtered)
filtered@int_colData@listData$reducedDims@listData$UMAP
vignette_cds@int_colData@listData$reducedDims@listData$UMAP
