#' Extract layer-0 embedding from Prov-GigaPath CSV
#'
#' @param pg_csv Path to a Prov-GigaPath CSV file
#' @param layer The embedding layer to extract (default 0)
#' @return Numeric vector of embeddings
#' @export
#' 
getEmbeddingLayer <- function(fnames, participant_ids = NULL, layer = NULL) {
  
  # If single file (character scalar), convert to vector for consistent handling
  single_file <- FALSE
  if (is.character(fnames) && length(fnames) == 1) {
    fnames <- list(fnames)
    single_file <- TRUE
    if (is.null(participant_ids)) participant_ids <- "sample1"
  }
  
  embeddings_list <- vector("list", length(fnames))
  
  # Use participant_ids as names if provided
  if (!is.null(participant_ids)) {
    names(embeddings_list) <- make.unique(participant_ids)
  } else {
    names(embeddings_list) <- paste0("sample", seq_along(fnames))
  }
  
  # Loop over files
  for (i in seq_along(fnames)) {
    if (file.exists(fnames[i])) {
      # Load CSV
      x <- readr::read_csv(fnames[[i]], show_col_types = FALSE)
      if (nrow(x) == 0) {
        embeddings_list[[i]] <- NA
        next
      }
      
      # Determine column index
      layer_ind <- if (is.null(layer)) {
        ncol(x)                  # last column
      } else if (is.character(layer)) {
        match(layer, colnames(x)) # find column by name
      } else {
        layer + 1                 # numeric index (0-based)
      }
      
      # Extract tensor string
      tensor_string <- x[[1, layer_ind]]
      # Clean string and convert to numeric
      clean_string <- gsub("tensor\\(\\[\\[|\\]\\]\\)", "", tensor_string)
      embeddings_list[[i]] <- as.numeric(unlist(strsplit(clean_string, ",\\s*")))
      
    } else {
      embeddings_list[[i]] <- NA
    }
  }
  
  # Remove fully NA entries
  embeddings_list <- embeddings_list[!sapply(embeddings_list, function(x) all(is.na(x)))]
  
  # If only one file, return the vector
  if (single_file) return(embeddings_list[[1]])
  
  # Convert list to matrix (pad with NA if needed)
  max_len <- max(sapply(embeddings_list, length))
  tensor_matrix <- t(sapply(embeddings_list, function(x) {
    if (length(x) < max_len) c(x, rep(NA, max_len - length(x))) else x
  }))
  
  rownames(tensor_matrix) <- names(embeddings_list)
  
  return(tensor_matrix)
}
