read_raw_PDUI =
  function (path){
    raw_PDUI = read.table(path, header = TRUE, stringsAsFactors = F)
    print(paste("number of genes in raw count matrix", nrow(raw_PDUI)))
    print(paste("number of cells in raw count matrix", ncol(raw_PDUI)))
    return(raw_PDUI)
  }

