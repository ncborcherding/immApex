#' Build Edit Distance Network Using Symmetric Deletion Lookup
#'
#' Constructs a weighted similarity network from biological sequences using a 
#' symmetric deletion lookup strategy combined with a banded edit-distance 
#' computation. The returned igraph object contains vertices representing
#' the input sequences and edges representing pairs of sequences whose 
#' edit distance is less than or equal to the specified threshold. The edge 
#' attribute \code{weight} stores the computed edit distance.
#'
#' This function supports both a character vector of sequences and a data frame
#' that contains at least three columns: \code{sequence}, \code{v.gene}, and
#' \code{j.gene}. Optionally, the function can filter candidate pairs by 
#' requiring matching \code{v.gene} and/or \code{j.gene} annotations 
#' (see \code{filter_v} and \code{filter_j}). 
#'
#' @param data A character vector of AIR sequences or a data frame with columns
#'   \code{sequence}, \code{v.gene}, and \code{j.gene}.
#' @param threshold An integer specifying the maximum allowed edit distance. Only
#'   pairs of sequences with an edit distance less than or equal to this value
#'   will be connected. Default is \code{2}.
#' @param filter.v Logical indicating whether to filter candidate pairs to only 
#' those that have matching \code{v.gene} family annotations. Default is 
#' \code{FALSE}.
#' @param filter.j Logical indicating whether to filter candidate pairs to only 
#' those that have matching \code{j.gene} family annotations. Default is 
#' \code{FALSE}.
#'
#' @return An igraph object representing the AIR similarity network. 
#' Vertices contain the original sequences (and gene annotations, if available),
#' and each edge has a \code{weight} attribute corresponding to the computed 
#' edit distance.
#'
#' @details
#' The function first calls a C++ routine (via Rcpp) to perform a symmetric 
#' deletion lookup,generating candidate pairs of sequences that might be within 
#' the specified edit distance.It then uses a banded dynamic programming 
#' algorithm (also implemented in C++) to compute the exact edit distance for 
#' each candidate pair. For sequences that are provided in a data frame, the 
#' function can optionally restrict comparisons to sequences with the same
#'  \code{v.gene} and/or \code{j.gene} values.
#'
#' @examples
#' sequences <- c("CASSLGTDTQYF", "CASSPGTDTQYF", "CASSLGNDTQYF", "CASRLGNDTQYF")
#' g <- buildNetwork(sequences, 
#'                   threshold = 2)
#' plot(g)
#'
#' df <- data.frame(
#'   sequence = c("CASSLGTDTQYF", "CASSPGTDTQYF", "CASSLGNDTQYF", "CASRLGNDTQYF"),
#'   v.gene = c("TRBV20", "TRBV20", "TRBV12", "TRBV20"),
#'   j.gene = c("TRBJ2-7", "TRBJ2-7", "TRBJ2-1", "TRBJ2-7")
#' )
#' g_df <- buildNetwork(df, 
#'                      threshold = 2, 
#'                      filter.v = TRUE, 
#'                      filter.j = TRUE)
#' plot(g_df)
#' 
#' @importFrom igraph graph_from_data_frame make_empty_graph V E
#' @importFrom Rcpp evalCpp
#' @export
buildNetwork <- function(data, 
                         threshold = 2, 
                         filter.v = FALSE, 
                         filter.j = FALSE) {
  
  #TODO V/J Gene Column Handling 
  #TODO V/J Gene Family Handling
  # Determine whether the input is a data frame or a character vector.
  if (is.data.frame(data)) {
    req_cols <- c("sequence", "v.gene", "j.gene")
    if (!all(req_cols %in% colnames(data))) {
      stop("Data frame must contain columns: sequence, v.gene, and j.gene")
    }
    sequences <- as.character(data$sequence)
    v_genes <- as.character(data[["v.gene"]])
    j_genes <- as.character(data[["j.gene"]])
  } else if (is.character(data)) {
    sequences <- data
    v_genes <- rep(NA, length(sequences))
    j_genes <- rep(NA, length(sequences))
  } else {
    stop("Input must be either a character vector or a data frame with the 
         required columns.")
  }
  
  n <- length(sequences)
  
  # Get candidate pairs using the symmetric deletion lookup (implemented in C++).
  candidate_pairs <- symmetric_deletion_lookup_cpp(sequences, threshold)
  
  # Post-filter candidates in parallel to verify edit distances and obtain edge weights.
  edge_df <- post_filter_candidates(candidate_pairs, 
                                    sequences, 
                                    v_genes, 
                                    j_genes,
                                    threshold, 
                                    filter_v, 
                                    filter_j)
  
  # Build vertices data frame to include all sequences.
  vertices_df <- data.frame(name = as.character(seq_len(n)), sequence = sequences,
                            stringsAsFactors = FALSE)
  if (!all(is.na(v_genes))) vertices_df$v.gene <- v_genes
  if (!all(is.na(j_genes))) vertices_df$j.gene <- j_genes
  
  # Create the igraph object with weighted edges.
  if (nrow(edge_df) == 0) {
    g <- make_empty_graph(n)
    V(g)$sequence <- sequences
    if (!all(is.na(v_genes))) V(g)$v.gene <- v_genes
    if (!all(is.na(j_genes))) V(g)$j.gene <- j_genes
  } else {
    g <- graph_from_data_frame(d = as.data.frame(edge_df), directed = FALSE, vertices = vertices_df)
    E(g)$weight <- as.numeric(edge_df$weight)
  }
  return(g)
}
