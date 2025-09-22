# -------------------------
# 1. GSEA plot function
# -------------------------
plot_gsea <- function(factor_id, gsea_collection = "hallmark", n_pathways = 8) {
  
  # Select GSEA results
  gsea_res <- switch(
    gsea_collection,
    "hallmark" = factor_interpretation[[factor_id]]$hallmark_top,
    "c2"       = factor_interpretation[[factor_id]]$c2_top,
    "c6"       = factor_interpretation[[factor_id]]$c6_top
  )
  
  if (is.null(gsea_res) || nrow(gsea_res) == 0) return(NULL)
  
  # Remove NAs
  gsea_res <- gsea_res[!is.na(gsea_res$pathway) & !is.na(gsea_res$NES) & !is.na(gsea_res$pval), ]
  if (nrow(gsea_res) == 0) return(NULL)
  
  # Select top pathways by p-value
  top_gsea <- head(gsea_res[order(gsea_res$pval), ], n_pathways)
  
  # Plot
  p <- ggplot(top_gsea, aes(x = NES, y = reorder(pathway, NES), fill = -log10(pval))) +
    geom_col() +
    scale_fill_gradient(low = "lightblue", high = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title = paste(factor_id, "GSEA:", gsea_collection),
         x = "Normalized Enrichment Score (NES)",
         y = "Pathway",
         fill = "-log10(pval)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}
