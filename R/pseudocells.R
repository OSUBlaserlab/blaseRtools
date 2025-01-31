bb_pseudocells <-
  function(tpm_matrix,
           n_cells,
           transcripts_per_cell,
           remove_genes = NULL) {
    exprs_list <- map(colnames(tpm_matrix),
                      .f = \(x,
                             dat = tpm_matrix,
                             remove = remove_genes) {
                        bulk <- x
                        vec <- dat[, bulk]
                        vec <- vec[!names(vec) %in% remove]
                        vec <- rep(names(vec), times = vec)

                        count_list <- map(.x = 1:n_cells, .f = \(x,
                                                                 cell = bulk,
                                                                 v = vec) {
                          cell <- paste0(cell, "_pc_", x)
                          v <-
                            sample(v,
                                   size = transcripts_per_cell,
                                   replace = FALSE)
                          tibble(original_rowname = v) |>
                            count(original_rowname, name = cell)
                        })

                        full_tbl <- reduce(count_list,
                                           .f = full_join,
                                           by = "original_rowname")
                        full_tbl[is.na(full_tbl)] <- 0

                        full_mtx <- as.matrix(full_tbl |>
                                                select(-original_rowname))
                        rownames(full_mtx) <-
                          full_tbl$original_rowname
                        as.sparse(full_mtx)

                      }) |> set_names(colnames(rounded_tpm_filtered))


    coldata_list <- map2(
      .x = exprs_list,
      .y = colnames(tpm_matrix),
      .f = \(x,y) {
        coldata <- tibble(cell_id = colnames(x),
                          dataset = y) |>
          as.data.frame()
        rownames(coldata) <- coldata$cell_id
        coldata$cell_id <- NULL
        coldata
      }
    )

    rowdata_list <- map(.x = exprs_list,
                        .f = \(x) {
                          data.frame(
                            id = rownames(x),
                            gene_short_name = rownames(x),
                            row.names = rownames(x)
                          )

                        })



    cell_line_cds_list <- pmap(
      .l = list(expr = exprs_list,
                col = coldata_list,
                row = rowdata_list),
      .f = \(expr, col, row) {
        new_cell_data_set(
          expression_data = expr,
          cell_metadata = col,
          gene_metadata = row
        )
      }
    )


    cell_line_cds <- combine_cds(cds_list = cell_line_cds_list)
    cell_line_cds <- preprocess_cds(cell_line_cds)
    cell_line_cds <- reduce_dimension(cell_line_cds)
    cell_line_cds

  }
