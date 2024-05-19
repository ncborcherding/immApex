percent.sequences <- function(input.sequences, 
                              max.length = NULL,
                              sequence.dictionary = amino.acids[1:20],
                              padding.symbol = ".") {
  sequence.dictionary <- c(sequence.dictionary, padding.symbol)
  if(is.null(max.length)) {
    max.length <- max(nchar(input.sequences))
  }
  
  #How to pad motifs
  print("Padding sequences...")
  padded_sequences <- .padded.strings(strings = input.sequences, 
                                      max.length = max.length,
                                      padded.token = padding.symbol,
                                      concatenate = TRUE)
  
  
  # Initialize the PSSM matrix with zeros
  pssm <- matrix(0, nrow = length(sequence.dictionary), ncol = max.length)
  rownames(pssm) <- sequence.dictionary
  
  # Convert sequences to a matrix
  sequence_matrix <- do.call(rbind, strsplit(unlist(padded_sequences), split = ""))
  
  # Count occurrences of each letter at each position using vectorized operations
  for (letter in sequence.dictionary) {
    pssm[letter, ] <- colSums(sequence_matrix == letter)
  }
  
  pssm <- pssm / length(input.sequences)
  
  colnames(pssm) <- paste0("Pos.", 1:max.length)
  
  return(pssm)
}
