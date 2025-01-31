break
# set up ref and qry ---------------------
data("hs_marrow_cds")
cds_ref <- hs_marrow_cds
load("data/cell_line_cds.rda")
cds_qry <- cell_line_cds

# Genes in reference.--------------------------
genes_ref <- row.names(cds_ref)
# Genes in query. -------------------------
genes_qry <- row.names(cds_qry)

# find a common name space-----------------------------
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
biomaRt::listAttributes(hsmart) |>
  as_tibble()


mapping <- getBM(
  attributes = c('ensembl_gene_id_version', 'ensembl_gene_id','entrezgene_id', 'hgnc_symbol'),
  mart = hsmart
) |>
  as_tibble() |>
  mutate(entrezgene_id = as.character(entrezgene_id))

qry_genes <- tibble(qry_genes = genes_qry) |>
  mutate(hgnc_symbol = str_remove(qry_genes, " .*"))

ref_genes <- tibble(ref_genes = genes_ref)

shared_namespace <- ref_genes |>
  mutate(ref_genes_short = str_remove(ref_genes, "\\..*")) |>
  inner_join(mapping, by = join_by(ref_genes_short == ensembl_gene_id)) |>
  filter(!is.na(hgnc_symbol)) |>
  left_join(qry_genes, by = join_by(hgnc_symbol)) |>
  filter(!is.na(qry_genes))

one_to_one_ref_genes <- shared_namespace |> count(ref_genes) |> filter(n ==1) |> pull(ref_genes)
one_to_one_qry_genes <- shared_namespace |> count(qry_genes) |> filter(n == 1) |> pull(qry_genes)

shared_namespace_one_to_one <- shared_namespace |>
  filter(ref_genes %in% one_to_one_ref_genes) |>
  filter(qry_genes %in% one_to_one_qry_genes)

cds_ref <- cds_ref[shared_namespace_one_to_one$ref_genes,]
cds_qry <- cds_qry[shared_namespace_one_to_one$qry_genes,]


# rebuild the cds objects with ensembl gene ids
cds_ref_exprs <- monocle3::exprs(cds_ref)
new_rownames_cds_ref_exprs <- tibble(old_rownames = rownames(cds_ref_exprs)) |>
  left_join(shared_namespace_one_to_one, by = join_by(old_rownames == ref_genes)) |>
  pull(ref_genes_short)

# sanity check
waldo::compare(str_remove(rownames(cds_ref_exprs), "\\..*"), new_rownames_cds_ref_exprs)

rownames(cds_ref_exprs) <- new_rownames_cds_ref_exprs

new_rowData_ref <- rowData(cds_ref)
rownames(new_rowData_ref) <- new_rownames_cds_ref_exprs
new_rowData_ref$id <- rownames(new_rowData_ref)

cds_ref_rebuilt <- new_cell_data_set(expression_data = cds_ref_exprs,
                                     cell_metadata = colData(cds_ref),
                                     gene_metadata = new_rowData_ref)

cds_qry_exprs <- monocle3::exprs(cds_qry)
new_rownames_cds_qry_exprs <- tibble(old_rownames = rownames(cds_qry_exprs)) |>
  left_join(shared_namespace_one_to_one, by = join_by(old_rownames == qry_genes)) |>
  pull(ref_genes_short)

rownames(cds_qry_exprs) <- new_rownames_cds_qry_exprs

new_rowData_qry <- rowData(cds_qry)
rownames(new_rowData_qry) <- new_rownames_cds_qry_exprs
new_rowData_qry
new_rowData_qry$id <- rownames(new_rowData_qry)


cds_qry_rebuilt <- new_cell_data_set(expression_data = cds_qry_exprs,
                                     cell_metadata = colData(cds_qry),
                                     gene_metadata = new_rowData_qry)



# Shared genes.
# Genes in rebuilt reference.--------------------------
genes_ref_rebuilt <- row.names(cds_ref_rebuilt)
# Genes in rebuilt query. -------------------------
genes_qry_rebuilt <- row.names(cds_qry_rebuilt)
genes_shared <- base::intersect(genes_ref_rebuilt, genes_qry_rebuilt)

# Remove non-shared genes and change the names --------------
cds_ref <- cds_ref_rebuilt[genes_shared,]
cds_qry <- cds_qry_rebuilt[genes_shared,]
cds_ref <- estimate_size_factors(cds_ref)
cds_qry <- estimate_size_factors(cds_qry)
cds_ref <- preprocess_cds(cds_ref, num_dim=100)
cds_ref <- reduce_dimension(cds_ref, build_nn_index=TRUE)
# Save the PCA and UMAP transform models for use with projection.
save_transform_models(cds_ref, 'cds_ref_test_models')
# Load the reference transform models into the query cds.
cds_qry <- load_transform_models(cds_qry, 'cds_ref_test_models')
# Apply the reference transform models to the query cds.
cds_qry <- preprocess_transform(cds_qry)
cds_qry <- reduce_dimension_transform(cds_qry)
plot_cells(cds_ref)
plot_cells(cds_qry)
# Label the data sets.
colData(cds_ref)[['data_set']] <- 'reference'
colData(cds_qry)[['data_set']] <- 'query'
# Combine the reference and query data sets.
cds_combined <- combine_cds(list(cds_ref, cds_qry),  keep_all_genes=TRUE, cell_names_unique=TRUE, keep_reduced_dims=TRUE)
plot_cells(cds_combined, color_cells_by='data_set')

bb_var_umap(cds_combined, "data_set")
bb_cellmeta(cds_combined) |> glimpse()

bb_var_umap(cds_combined, "celltype.l2", value_to_highlight = c("HSC", "LMPP", "GMP"))
bb_var_umap(cds_combined, "celltype.l2")
bb_var_umap(cds_combined, "celltype.l2", overwrite_labels = TRUE) +
  bb_var_umap(cds_combined, "cell_line")

cds_qry_lab_xfr <- transfer_cell_labels(cds_qry, reduction_method='UMAP', ref_coldata=colData(cds_ref), ref_column_name='celltype.l2', query_column_name='cell_type_xfr', transform_models_dir='cds_ref_test_models')
cds_qry_lab_fix <- fix_missing_cell_labels(cds_qry_lab_xfr, reduction_method='UMAP', from_column_name='cell_type_xfr', to_column_name='cell_type_fix')

bb_var_umap(cds_qry_lab_fix, "cell_type_fix")

cds_combined <- bb_cellmeta(cds_ref) |>
  select(cell_id, merged_celltype.l2 = celltype.l2) |>
  bind_rows(bb_cellmeta(cds_qry_lab_fix) |>
              select(cell_id, merged_celltype.l2 = cell_type_fix)) |>
  bb_tbl_to_coldata(cds_combined, min_tbl = _)

bb_var_umap(cds_combined, "merged_celltype.l2", facet_by = "data_set", overwrite_labels = TRUE)

bm_cell_line_cds <- cds_combined


library(blaseRdata)
vignette_cds |> bb_cellmeta()
test <- bb_monocle_(cds_qry = vignette_cds,
            cds_ref = vignette_cds,
            labels = "seurat_celltype_l2",
            use_aligned = FALSE,
            suffix = "_ref")
test <- SingleCellExperiment::reducedDim(vignette_cds)
test
waldo::compare(test, SingleCellExperiment::reducedDims(vignette_cds)$PCA)

test_cds <- vignette_cds
reducedDim(test_cds) <- reducedDims(vignette_cds)$Aligned
reducedDims(test_cds)
waldo::compare(SingleCellExperiment::reducedDims(test_cds)$PCA, SingleCellExperiment::reducedDims(test_cds)$Aligned)
waldo::compare(SingleCellExperiment::reducedDim(vignette_cds)$PCA, SingleCellExperiment::reducedDims(vignette_cds)$Aligned)

test <- bb_monocle_(
  cds_qry = vignette_cds,
  cds_ref = vignette_cds,
  labels = "seurat_celltype_l2",
  use_aligned = FALSE,
  suffix = "_ref"
)

bb_cellmeta(vignette_cds)

vignette_cds_head <- filter_cds(vignette_cds, cells = bb_cellmeta(vignette_cds) |> slice_head(n = 500))
vignette_cds_tail <- filter_cds(vignette_cds, cells = bb_cellmeta(vignette_cds) |> slice_tail(n = 500))
test1 <- bb_monocle_project(
  cds_qry = vignette_cds_head,
  cds_ref = vignette_cds_tail,
  labels = c("seurat_celltype_l2", "seurat_celltype_l1"),
  use_aligned = TRUE,
  suffix = "_ref"
)

test2 <- bb_monocle_anno(
  cds_qry = vignette_cds_head,
  cds_ref = vignette_cds_tail,
  labels = c("seurat_celltype_l2", "seurat_celltype_l1"),
  use_aligned = TRUE,
  suffix = "_ref"
)
bb_cellmeta(test1) |> glimpse()
bb_var_umap(test1, "data_set")
bb_var_umap(test1, "merged_seurat_celltype_l2")
bb_var_umap(test1, "merged_seurat_celltype_l1")

bb_cellmeta(test2) |> glimpse()

a <- "appple"
b <- "banana"
c <- "carrot"

c(a, b, c)
c("a", "b", "c")
!!sym("a")
dplyr::sy

debug_tbl <- tibble(
  cell_id = rep("something", times = 10),
  a = rep("apple", times = 10),
       b = rep("banana", times = 10),
       c = rep("carrot", times = 10))



bb_monocle_debug <-
  function(new_cellmeta = debug_tbl |> mutate(a_ref = b),
           labels,
           suffix = "_ref",
           use_aligned = TRUE) {




    # merge the transferred colData columns
    new_cellmeta <- map(labels, .f = \(x,
                                       dat = new_cellmeta,
                                       suf = suffix) {
      dat |>
        select(cell_id, !!sym(paste0("merged_", x)) := all_of(paste0(x, suffix)))
    }) |> reduce(.f = left_join)
    new_cellmeta
    # ref_cellmeta <- bb_cellmeta(cds_ref) |>
    #   dplyr::select(cell_id, tidyr::matches(labels))
    # ref_cellmeta[[paste0("merged_", labels)]] <-
    #   ref_cellmeta[[labels]]
    # new_cellmeta <- dplyr::bind_rows(new_cellmeta, ref_cellmeta)
    # new_cellmeta <-
    #   new_cellmeta[c("cell_id", paste0("merged_", labels))]
    #
    # cds_combined <-
    #   bb_tbl_to_coldata(cds_combined, min_tbl = new_cellmeta)


  }
bb_monocle_debug(labels = c("a", "b"))
bb_monocle_debug(labels = c("a"))

x <- "barcode"
n <- 1
while (x %in% colnames(colData(vignette_cds))) {
  x <- paste0(x, ".", n)
  n <- n + 1
}
x

test_cond <- c("barcode", "barcode.1", "barcode.2", "barcode.3")
x <- "barcode"
n <- 1
while (x %in% test_cond) {
  x <- paste0(str_remove(x,"\\..*"), ".", n)
  n <- n + 1
}
x


purrr::map2(.x = c("A", "B", "C"),
            .y = c("1"),
            .f = \(x,y) {
              paste0(x, y)
            }
            )
