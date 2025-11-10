#' Run GSVA/ssGSEA on a selected MOFA view
#'
#' @param mofa_object Trained MOFA object
#' @param pathways Named list of gene sets (e.g., from get_gene_sets())
#' @param view Optional. Name of the MOFA view to use. If NULL, shows available views and prompts the user
#' @param method GSVA method (default: "ssgsea")
#' @param verbose Logical, passed to gsva
#'
#' @return Matrix of enrichment scores (pathways x samples)
#' @export
run_gsva_mofa <- function(mofa_object, pathways, view = NULL, verbose = TRUE) {
  # List available views
  available_views <- names(mofa_object@data)
  
  # Prompt user if view not provided
  if (is.null(view)) {
    message("Available views: ", paste(available_views, collapse = ", "))
    view <- readline(prompt = "Enter view to use: ")
  }
  if (!view %in% available_views) {
    stop("View '", view, "' not found. Available views: ", paste(available_views, collapse = ", "))
  }
  
  # Extract expression matrix
  expr_list <- get_data(mofa_object, view = view)[[1]]
  expr_mat <- as.matrix(expr_list[[1]])  # assumes single group per view
  
  # Create GSVA parameters
  ssgsea_params <- ssgseaParam(expr_mat, pathways)
  
  # Run GSVA
  gsva_scores <- gsva(ssgsea_params, verbose = verbose)
  
  return(gsva_scores)
}