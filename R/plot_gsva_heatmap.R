# -------------------------
# GSVA heatmap
# -------------------------
plot_gsva_heatmap <- function(gsva_collection = "hallmark", top_n_pathways = 20) {
  
  gsva_list <- switch(
    gsva_collection,
    "hallmark" = hallmark_interpret,
    "c2"       = c2_interpret,
    "c6"       = c6_interpret
  )
  
  gsva_mat <- sapply(names(gsva_list), function(factor_id) {
    corr <- gsva_list[[factor_id]]$positive_cor
    names(corr) <- gsva_list[[factor_id]]$positive_pathways
    corr <- corr[!is.na(names(corr)) & !is.na(corr)]
    head_corr <- head(corr, top_n_pathways)
    if(length(head_corr) < top_n_pathways) head_corr <- c(head_corr, rep(NA, top_n_pathways - length(head_corr)))
    names(head_corr) <- names(corr)[1:length(head_corr)]
    return(head_corr)
  })
  
  gsva_mat <- as.matrix(gsva_mat)
  
  # Remove rows/cols all NA
  gsva_mat <- gsva_mat[rowSums(!is.na(gsva_mat)) > 0, , drop = FALSE]
  gsva_mat <- gsva_mat[, colSums(!is.na(gsva_mat)) > 0, drop = FALSE]
  
  pheatmap(gsva_mat,
           color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
           main = paste("GSVA Correlation Heatmap -", gsva_collection),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 8,
           fontsize_col = 10,
           na_col = "grey")
}