#' Retrieve MSigDB Gene Sets
#'
#' Fetches MSigDB collections (Hallmark, C2, C6) and returns them as lists
#' compatible with `fgsea` or `GSVA`.
#'
#' @param species Character. Species name (default: "Homo sapiens").
#' @param collections Character vector. Any of "H", "C2", "C6".
#' @param verbose Logical. If TRUE, prints messages.
#'
#' @return Named list of gene set collections. Each element is a list of pathways.
#' @examples
#' gene_sets <- get_gene_sets(collections = c("H", "C2", "C6"))
#' names(gene_sets)
#' length(gene_sets$hallmark)
#' @export
get_gene_sets <- function(species = "Homo sapiens",
                          collections = c("H", "C2", "C6"),
                          verbose = TRUE) {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required but not installed.")
  
  pathways_list <- list()
  
  for (col in collections) {
    if (verbose) message("Fetching ", col, " gene sets...")
    
    msig <- msigdbr::msigdbr(species = species, collection = col)
    
    # ensure both gene_symbol and gs_name exist
    if (!all(c("gene_symbol", "gs_name") %in% colnames(msig))) {
      stop("msigdbr output missing expected columns.")
    }
    
    gs_split <- split(x = msig$gene_symbol, f = msig$gs_name)
    
    # assign human-readable name
    if (col == "H") {
      pathways_list$hallmark <- gs_split
    } else if (col == "C2") {
      pathways_list$c2 <- gs_split
    } else if (col == "C6") {
      pathways_list$c6 <- gs_split
    } else {
      warning("Unknown collection: ", col, ", skipping.")
    }
  }
  
  if (verbose) message("Pathway retrieval complete.")
  return(pathways_list)
}