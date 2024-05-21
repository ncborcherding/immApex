substitution.matrix <- function(input.sequences = NULL, 
                                normalize = FALSE,
                                sequence.dictionary = amino.acids[1:20]) {
  # Initialize the substitution matrix with zeros
  n <- length(sequence.dictionary)
  substitution_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(substitution_matrix) <- sequence.dictionary
  colnames(substitution_matrix) <- sequence.dictionary
  
  # Create a lookup for dictionary positions
  dict_lookup <- setNames(seq_along(sequence.dictionary), sequence.dictionary)
  
  # Iterate through each pair of sequences
  for (i in 1:length(input.sequences)) {
    seq1 <- strsplit(input.sequences[i], "")[[1]]
    len1 <- length(seq1)
    for (j in (i + 1):length(input.sequences)) {
      if (j > length(input.sequences)) next # Check to ensure j is within bounds
      
      seq2 <- strsplit(input.sequences[j], "")[[1]]
      len2 <- length(seq2)
      min_length <- min(len1, len2)
      
      # Count substitutions
      for (pos in 1:min_length) {
        letter1 <- seq1[pos]
        letter2 <- seq2[pos]
        
        # Check if both letters are in the dictionary
        if (!is.null(dict_lookup[[letter1]]) && !is.null(dict_lookup[[letter2]])) {
          idx1 <- dict_lookup[[letter1]]
          idx2 <- dict_lookup[[letter2]]
          substitution_matrix[idx1, idx2] <- substitution_matrix[idx1, idx2] + 1
          substitution_matrix[idx2, idx1] <- substitution_matrix[idx2, idx1] + 1
        }
      }
    }
  }
  
  
  # Normalize the substitution matrix if required
  if (normalize) {
    substitution_matrix <- substitution_matrix / sum(substitution_matrix)
  }
  
  return(substitution_matrix)
}