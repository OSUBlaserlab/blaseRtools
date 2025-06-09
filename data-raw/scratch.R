devtools::load_all()


filter_cds(vignette_cds, cells = bb_cellmeta(vignette_cds) |> dplyr::filter(sample == "chromium_X"))
filter_cds(vignette_cds, cells = bb_cellmeta(vignette_cds) |> dplyr::filter(sample == "chromium_Y"))
filter_cds(vignette_cds, cells = "all")

filter_cds(vignette_cds, genes = bb_rowmeta(vignette_cds) |> dplyr::filter(gene_short_name == "CD8A"))
filter_cds(vignette_cds, genes = bb_rowmeta(vignette_cds) |> dplyr::filter(gene_short_name == "CD8Z"))
filter_cds(vignette_cds, genes = "all")
