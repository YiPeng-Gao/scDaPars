Processing_raw_PDUI_matrix =
  function(raw_PDUI, filter_gene_thre, filter_cell_thre)
  {
    raw_PDUI_matrix_sc = raw_PDUI[,5:ncol(raw_PDUI)]
    genenames = sapply(strsplit(raw_PDUI$Gene, "[|]"), `[`, 1)
    row.names(raw_PDUI_matrix_sc) = genenames
    raw_PDUI_matrix_sc = raw_PDUI_matrix_sc[order(row.names(raw_PDUI_matrix_sc)),]
    raw_PDUI_matrix_sc = raw_PDUI_matrix_sc[which(rowSums(!is.na(raw_PDUI_matrix_sc)) >= ncol(raw_PDUI_matrix_sc) * filter_gene_thre),]
    pass_cells_Dapars= names(which(colSums(!is.na(raw_PDUI_matrix_sc)) >= nrow(raw_PDUI_matrix_sc) * filter_cell_thre))
    raw_PDUI_matrix_sc = subset(raw_PDUI_matrix_sc, select = which(colnames(raw_PDUI_matrix_sc) %in% pass_cells_Dapars))
    raw_PDUI_matrix_sNA = raw_PDUI_matrix_sc
    raw_PDUI_matrix_sNA[is.na(raw_PDUI_matrix_sNA)] = 0
    return(list(raw_PDUI_matrix_sc = raw_PDUI_matrix_sc, raw_PDUI_matrix_sNA = raw_PDUI_matrix_sNA))
  }
