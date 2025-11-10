#' Run FGSEA for all MOFA factors
#'
#' @description
#' Performs FGSEA enrichment across all MOFA factors for a given list of gene sets.
#'
#' @param mofa_object Trained MOFA object.
#' @param pathways Named list of gene sets. Each element is a character vector of gene symbols.
#' @param view Character. MOFA view to use (default = "RNASeq").
#' @param minSize Minimum pathway size (default 15).
#' @param maxSize Maximum pathway size (default 500).
#' @param top_n Integer. Number of top pathways to return per factor (default 20).
#'
#' @return A named list of FGSEA results per factor.
#' @examples
#' \dontrun{
#' fgsea_results <- run_fgsea_all_factors(
#'   mofa_object.trained,
#'   pathways = gene_sets$hallmark
#' )
#' }
#' @export
run_fgsea_all_factors <- function(mofa_object,
                                  pathways,
                                  view = "RNASeq",
                                  minSize = 15,
                                  maxSize = 500,
                                  top_n = 20) {
  stopifnot(requireNamespace("fgsea", quietly = TRUE))
  
  if (!is.list(pathways)) {
    stop("pathways must be a named list of gene sets")
  }
  
  weights_df <- get_weights(mofa_object, views = view, factors = "all", as.data.frame = TRUE)
  res_list <- list()
  
  for (factor_id in unique(weights_df$factor)) {
    fw <- subset(weights_df, factor == factor_id)
    ranks <- fw$value
    names(ranks) <- fw$feature
    ranks <- sort(ranks, decreasing = TRUE)
    
    fg <- fgsea::fgsea(pathways = pathways, stats = ranks, minSize = minSize, maxSize = maxSize)
    fg <- fg[order(fg$pval), ]
    res_list[[factor_id]] <- head(fg, top_n)
  }
  
  return(res_list)
}