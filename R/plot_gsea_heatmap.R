# -------------------------
# GSEA heatmap
# -------------------------
plot_gsea_heatmap <- function(gsea_collection = "hallmark", top_n_pathways = 20) {
  
  gsea_mat <- sapply(names(factor_interpretation), function(factor_id) {
    
    gsea_res <- switch(
      gsea_collection,
      "hallmark" = factor_interpretation[[factor_id]]$hallmark_top,
      "c2"       = factor_interpretation[[factor_id]]$c2_top,
      "c6"       = factor_interpretation[[factor_id]]$c6_top
    )
    
    if (is.null(gsea_res) || nrow(gsea_res) == 0) return(rep(NA, top_n_pathways))
    
    gsea_res <- gsea_res[!is.na(gsea_res$pathway) & !is.na(gsea_res$NES), ]
    gsea_res <- gsea_res[order(gsea_res$pval), ]
    
    top_res <- head(gsea_res, top_n_pathways)
    NES_vec <- top_res$NES
    names(NES_vec) <- top_res$pathway
    if(length(NES_vec) < top_n_pathways) NES_vec <- c(NES_vec, rep(NA, top_n_pathways - length(NES_vec)))
    
    return(NES_vec)
  })
  
  gsea_mat <- as.matrix(gsea_mat)
  
  # Remove rows (pathways) that are all NA
  gsea_mat <- gsea_mat[rowSums(!is.na(gsea_mat)) > 0, , drop = FALSE]
  gsea_mat <- gsea_mat[, colSums(!is.na(gsea_mat)) > 0, drop = FALSE]
  
  pheatmap(gsea_mat,
           color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
           main = paste("GSEA NES Heatmap -", gsea_collection),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 8,
           fontsize_col = 10,
           na_col = "grey")
}
