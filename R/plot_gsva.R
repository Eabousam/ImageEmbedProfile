
# -------------------------
# 2. GSVA plot function
# -------------------------
plot_gsva <- function(factor_id, gsva_collection = "hallmark", n_pathways = 8) {
  
  # Select GSVA interpretation
  gsva_list <- switch(
    gsva_collection,
    "hallmark" = hallmark_interpret[[factor_id]],
    "c2"       = c2_interpret[[factor_id]],
    "c6"       = c6_interpret[[factor_id]]
  )
  
  if (is.null(gsva_list)) return(NULL)
  
  # Remove NAs
  pos_idx <- which(!is.na(gsva_list$positive_pathways) & !is.na(gsva_list$positive_cor))
  neg_idx <- which(!is.na(gsva_list$negative_pathways) & !is.na(gsva_list$negative_cor))
  
  top_pos <- head(gsva_list$positive_pathways[pos_idx], n_pathways)
  top_neg <- head(gsva_list$negative_pathways[neg_idx], n_pathways)
  cor_pos <- head(gsva_list$positive_cor[pos_idx], n_pathways)
  cor_neg <- head(gsva_list$negative_cor[neg_idx], n_pathways)
  
  # Skip if empty
  if (length(top_pos) + length(top_neg) == 0) return(NULL)
  
  # Combine into data.frame
  plot_data <- data.frame(
    pathway = c(top_pos, top_neg),
    cor = c(cor_pos, cor_neg),
    type = c(rep("Positive", length(top_pos)), rep("Negative", length(top_neg)))
  )
  
  # Plot
  p <- ggplot(plot_data, aes(x = cor, y = reorder(pathway, cor), fill = type)) +
    geom_col() +
    scale_fill_manual(values = c("Positive" = "red", "Negative" = "darkgreen")) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title = paste(factor_id, "GSVA:", gsva_collection),
         x = "Spearman correlation",
         y = "Pathway",
         fill = "Correlation") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}