#' Position Probability Matrix for Amino Acid or Nucleotide Sequences
#' 
#' Use this to generate a position-probability or weight matrix 
#' for a set of given sequences. 
#' 
#' @examples
#' new.sequences <- generate.sequences(prefix.motif = "CAS",
#'                                     suffix.motif = "YF",
#'                                     number.of.sequences = 100,
#'                                     min.length = 8,
#'                                     max.length = 16)
#'                           
#' PPM.matrix <- probability.matrix(new.sequences)
#'                         
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param max.length Additional length to pad, NULL will pad sequences 
#' to the max length of input.sequences
#' @param convert.PWM Convert the matrix into a positional weight matrix 
#' using log likelihood
#' @param background.frequencies Provide amino acid or nucleotide frequencies
#' for the positional weight matrix. If NULL, assumes uniform likelihood.
#' @param sequence.dictionary The letters to use in sequence generation 
#' (default are all amino acids)
#' @param padding.symbol Symbol to use for padding at the end of sequences
#' 
#' @importFrom stats median
#' 
#' @export
#' @return A matrix with position specific probabilities or weights
probability.matrix <- function(input.sequences, 
                               max.length = NULL,
                               convert.PWM = FALSE,
                               background.frequencies = NULL,
                               sequence.dictionary = amino.acids[1:20],
                               padding.symbol = ".") {
  sequence.dictionary <- c(sequence.dictionary, padding.symbol)
  if(!is.null(background.frequencies)) {
    if(length(background.frequencies) != length(sequence.dictionary)-1) {
      stop("Please ensure the background.frequencies match the length of the sequence.dictionary")
    }
    #Adding a frequency for padded values
    background.frequencies <- c(background.frequencies, median(background.frequencies))
    names(background.frequencies) <- sequence.dictionary
  }
  
  
  if(is.null(max.length)) {
    max.length <- max(nchar(input.sequences))
  }
  
  print("Padding sequences...")
  padded_sequences <- .padded.strings(strings = input.sequences, 
                                      max.length = max.length,
                                      padded.token = padding.symbol,
                                      concatenate = TRUE)
  
  print("Calculating Positional Probabilites for sequences...")
  # Initialize the PSSM matrix with zeros
  position_matrix <- matrix(0, nrow = length(sequence.dictionary), ncol = max.length)
  rownames(position_matrix) <- sequence.dictionary
  
  # Convert sequences to a matrix
  sequence_matrix <- do.call(rbind, strsplit(unlist(padded_sequences), split = ""))
  
  # Count occurrences of each letter at each position using vectorized operations
  for (letter in sequence.dictionary) {
    position_matrix[letter, ] <- colSums(sequence_matrix == letter)
  }
  
  #Normalizing
  if(convert.PWM) {
    position_matrix <- position_matrix + 1
    position_matrix<- position_matrix / (length(input.sequences) + length(sequence.dictionary) * 1)
  } else {
    position_matrix <- position_matrix / length(input.sequences)
  }
  
  colnames(position_matrix) <- paste0("Pos.", seq_len(max.length))
  
  if(convert.PWM) {
    print("Converting to Liklihoods for a PWM...")
    # Calculate log-likelihood PSSM
    if (is.null(background.frequencies)) {
      # If no background frequencies provided, assume uniform distribution
      background.frequencies <- rep(1 / length(sequence.dictionary), length(sequence.dictionary))
      names(background.frequencies) <- sequence.dictionary
    } 
    position_matrix <- log2(position_matrix / background.frequencies[rownames(position_matrix)])
  }
  
  return(position_matrix)
}
