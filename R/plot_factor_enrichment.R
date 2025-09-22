plot_factor_enrichment <- function(factor_id, 
                                   pathway_collection = "hallmark", 
                                   top_n = 1,
                                   all_weights) {
  # Extract weights for the factor
  factor_weights <- all_weights[all_weights$factor == factor_id, ]
  if(nrow(factor_weights) == 0) {
    message("No weights found for factor: ", factor_id)
    return(NULL)
  }
  
  # Create ranked vector
  ranks <- sort(setNames(as.numeric(factor_weights$value), as.character(factor_weights$feature)),
                decreasing = TRUE)
  
  # Select pathways collection
  pathways <- switch(pathway_collection,
                     "hallmark" = pathways_hallmark,
                     "c2"       = pathways_c2,
                     "c6"       = pathways_c6,
                     stop("Invalid pathway collection"))
  
  # Filter pathways to those overlapping the ranked genes
  valid_pathways <- pathways[sapply(pathways, function(p) any(p %in% names(ranks)))]
  if(length(valid_pathways) == 0) {
    message("No overlapping pathways found for factor: ", factor_id)
    return(NULL)
  }
  
  # Run FGSEA
  fgsea_res <- fgsea(pathways = valid_pathways, stats = ranks, minSize = 15, maxSize = 500)
  fgsea_res <- fgsea_res[!is.na(fgsea_res$pathway), ]
  if(nrow(fgsea_res) == 0) {
    message("FGSEA returned no results for factor: ", factor_id)
    return(NULL)
  }
  
  # Select top pathways by p-value
  fgsea_top <- head(fgsea_res[order(fgsea_res$pval), ], top_n)
  if(nrow(fgsea_top) == 0) {
    message("No top pathways available for plotting")
    return(NULL)
  }
  
  # Plot enrichment for top pathway(s)
  enrichment_plots <- lapply(fgsea_top$pathway, function(path) {
    plotEnrichment(valid_pathways[[path]], ranks) +
      labs(title = paste("Factor", factor_id, "-", path))
  })
  
  return(enrichment_plots)
}

