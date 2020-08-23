Find_Neighboring_Cells =
  function(raw_PDUI_sNA)
  {
    var_thre = 0.4 # only keep PCs that can explain at least 40% of variance
    pca = prcomp(t(raw_PDUI_sNA))
    eigs = (pca$sdev)^2
    var_cum = cumsum(eigs)/sum(eigs)
    npc = which.max(var_cum > var_thre)
    mat_pcs = t(pca$x[,1:npc]) # columns are cells

    # detect outliers cells
    dist_cell = as.data.frame(as.matrix(dist(t(mat_pcs), diag = T, upper = T)))
    min_dist = sapply(1:ncol(mat_pcs), function(i){min(dist_cell[i, -i])}) # find minimum distance
    iqr = quantile(min_dist, 0.75) - quantile(min_dist, 0.25)
    outliers = which(min_dist > 1.5 * iqr + quantile(min_dist, 0.75))
    non_out = setdiff(1:ncol(mat_pcs), outliers)
    outliers = colnames(raw_PDUI_sNA)[outliers]
    mat_pcs = mat_pcs[,non_out]

    # Graph based clustering to find potential simialr cells (without specify # of clusters)
    knn.info = RANN::nn2(t(mat_pcs), k=30)
    knn = knn.info$nn.idx
    adj = matrix(0, ncol(mat_pcs), ncol(mat_pcs))
    rownames(adj) = colnames(adj) = colnames(mat_pcs)
    for (i in seq_len(ncol(mat_pcs))) {
      adj[i,colnames(mat_pcs)[knn[i,]]] = 1
    }
    g <- igraph::graph.adjacency(adj, mode="undirected")
    g <- simplify(g)
    km = igraph::cluster_walktrap(g)
    com = km$membership
    names(com) = km$names
    num_clusters = max(km$membership)
    cluster_info = data.frame(km$membership)
    row.names(cluster_info) = km$names
    return(list(cluster_info = cluster_info, num_clusters = num_clusters, outliers = outliers))
  }
