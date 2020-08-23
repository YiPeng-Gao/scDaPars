impute_nnls =
  function(ngenes, cellid, DaPars_table_raw, DaPars_table_softImpute, geneid_drop,
                       geneid_obs, nbs, distc)
  {
    yobs = DaPars_table_raw[ ,cellid]
    if (length(geneid_drop) == 0 | length(geneid_drop) == ngenes) {
      return(yobs) }

    Yimpute = rep(0, ngenes)
    XNrobust = as.matrix(DaPars_table_raw[geneid_obs, nbs])
    Yrobust = as.matrix(DaPars_table_raw[geneid_obs, cellid])
    XNimpute = as.matrix(DaPars_table_softImpute[geneid_drop, nbs])

    num_thre = 500
    if(ncol(XNrobust) >= min(num_thre, nrow(XNrobust))){
      if (num_thre >= nrow(XNrobust)){
        new_thre = round((2*nrow(XNrobust)/3))
      }else{ new_thre = num_thre}
      filterid = order(distc[cellid, -cellid])[1: new_thre]
      XNrobust = XNrobust[, filterid, drop = FALSE]
      XNimpute = XNimpute[, filterid, drop = FALSE]
    }

    set.seed(cellid)
    nnls = penalized(Yrobust, penalized = XNrobust, unpenalized = ~0,
                     positive = TRUE, lambda1 = 0, lambda2 = 0,
                     maxiter = 3000, trace = FALSE)
    Ynew = penalized::predict(nnls, penalized = XNimpute, unpenalized = ~0)[,1]

    Yimpute[geneid_drop] = Ynew
    Yimpute[geneid_obs] = yobs[geneid_obs]
    maxobs = rep(1,ngenes)
    Yimpute[Yimpute > maxobs] = maxobs[Yimpute > maxobs]

    return(Yimpute)
  }


Impuation_model =
  function(raw_PDUI_sNA, raw_PDUI_sc, num_clusters, cluster_info)
  {
    PDUI_sc_imputed = data.frame(row.names = row.names(raw_PDUI_sc))
    for(i in 1:num_clusters){
      cluster_cells = row.names(cluster_info)[which(cluster_info$km.membership == i)]
      cluster_raw_PDUI_sc = raw_PDUI_sc[, which(colnames(raw_PDUI_sc) %in% cluster_cells)]
      cluster_raw_PDUI_sNA = raw_PDUI_sNA[, which(colnames(raw_PDUI_sNA) %in% cluster_cells)]

      #identify robust and dropout genes
      setA = lapply(1:ncol(cluster_raw_PDUI_sc), function(cell) which(is.na(cluster_raw_PDUI_sc[,cell]))) #dropouts
      setB = lapply(1:ncol(cluster_raw_PDUI_sc), function(cell) which(!is.na(cluster_raw_PDUI_sc[,cell]))) #robust genes

      #calculae DaPars cell distance
      dist_cells_list = lapply(1:ncol(cluster_raw_PDUI_sNA), function(id1){
        d = sapply(1:id1, function(id2){
          sse = sum((cluster_raw_PDUI_sNA[, id1] - cluster_raw_PDUI_sNA[, id2])^2, na.rm = T)
          sqrt(sse)
        })
        return(c(d, rep(0, ncol(cluster_raw_PDUI_sNA)-id1)))
      })
      cell_dist = matrix(0, nrow = ncol(cluster_raw_PDUI_sNA), ncol = ncol(cluster_raw_PDUI_sNA))
      for(cellid in 1:ncol(cluster_raw_PDUI_sNA)){cell_dist[cellid, ] = dist_cells_list[[cellid]]}
      cell_dist = cell_dist + t(cell_dist)

      cluster_imputation.res = foreach(cellid = 1:ncol(cluster_raw_PDUI_sc), .packages = c("penalized"),
                                       .combine = cbind, .export = c("impute_nnls")) %do% {
                                         if (cellid %% 10 == 0) {gc()}
                                         if (cellid %% 100 == 0) {print(cellid)}
                                         nbs = setdiff(1:ncol(cluster_raw_PDUI_sc), cellid)
                                         if (length(nbs) == 0) {return(NULL)}
                                         geneid_drop = setA[[cellid]]
                                         geneid_obs = setB[[cellid]]
                                         y = try(impute_nnls(nrow(cluster_raw_PDUI_sc), cellid, cluster_raw_PDUI_sNA, cluster_raw_PDUI_sNA,
                                                             geneid_drop, geneid_obs, nbs, distc = cell_dist),
                                                 silent = TRUE)
                                         if (class(y) == "try-error") {
                                           y = cluster_raw_PDUI_sNA[, cellid, drop = FALSE]
                                         }
                                         return(y)
                                       }
      rownames(cluster_imputation.res) = rownames(cluster_raw_PDUI_sc)
      colnames(cluster_imputation.res) = colnames(cluster_raw_PDUI_sc)

      PDUI_sc_imputed = cbind(PDUI_sc_imputed, cluster_imputation.res)
    }
    return(PDUI_sc_imputed)
  }
