#' Adjacency matrix from amino acid or nucleotide sequences
#' 
#' Calculate frequency of adjacency between residues
#' along a set of biological sequences.
#' 
#' @examples
#' new.sequences <- generateSequences(prefix.motif = "CAS",
#'                                    suffix.motif = "YF",
#'                                    number.of.sequences = 100,
#'                                    min.length = 8,
#'                                    max.length = 16)
#'                           
#' adj.matrix <- adjacencyMatrix(new.sequences,
#'                               normalize = TRUE)
#'                         
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param normalize Return the values as a function of total number of 
#' residues (\strong{TRUE}) or frequencies (\strong{FALSE})
#' @param sequence.dictionary The letters to use in sequence generation 
#' (default are all amino acids)
#' 
#' @export
#' @return Adjacency matrix based on input.sequences.
adjacencyMatrix <- function(input.sequences = NULL, 
                            normalize = TRUE,
                            sequence.dictionary = amino.acids) {
  
  # Initialize the adjacency matrix with zeros
  n <- length(sequence.dictionary)
  adjacency_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(adjacency_matrix) <- sequence.dictionary
  colnames(adjacency_matrix) <- sequence.dictionary
  
  # Create a lookup for dictionary positions
  dict_lookup <- setNames(seq_along(sequence.dictionary), sequence.dictionary)
  
  # Iterate through each sequence
  for (seq in input.sequences) {
    seq_chars <- strsplit(seq, "")[[1]]
    len <- length(seq_chars)
    
    # Count adjacency
    for (pos in seq_len(len - 1)) {
      letter1 <- seq_chars[pos]
      letter2 <- seq_chars[pos + 1]
      
      # Check if both letters are in the dictionary
      if (!is.null(dict_lookup[[letter1]]) && !is.null(dict_lookup[[letter2]])) {
        idx1 <- dict_lookup[[letter1]]
        idx2 <- dict_lookup[[letter2]]
        adjacency_matrix[idx1, idx2] <- adjacency_matrix[idx1, idx2] + 1
        adjacency_matrix[idx2, idx1] <- adjacency_matrix[idx2, idx1] + 1
      }
    }
  }
  
  # Normalize the adjacency matrix if required
  if (normalize) {
    adjacency_matrix <- adjacency_matrix / sum(adjacency_matrix)
  }
  
  return(adjacency_matrix)
}