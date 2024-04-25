#' One Hot Encoder from Amino Acid String
#' 
#' Use this to transform amino acid sequences a
#' one hot encoding of the sequence.
#' 
#' @examples
#' new.sequences <- generate.sequences(prefix.motif = "CAS",
#'                                     suffix.motif = "YF",
#'                                     number.of.sequences = 100,
#'                                     min.length = 8,
#'                                     max.length = 16)
#'                           
#'sequence.matrix <- one.hot.encoder(new.sequences, 
#'                                   convert.to.matrix = TRUE)
#'                         
#' @param input.sequences The amino acid sequences to use
#' @param max.length Additional length to pad, NULL will pad sequences 
#' to the max length of input.sequences
#' @param convert.to.matrix Return a matrix (TRUE) or a 3D array (FALSE)
#' 
#' @export
#' @return Tokenized amino acid sequences in a matrix of vector
one.hot.encoder <- function(input.sequences, 
                            max.length = NULL,
                            convert.to.matrix = TRUE) {
  
  char_set <- c(amino.acids, ".")
  # Create a mapping of amino acids to integers
  char_to_int <- setNames(seq_along(char_set), char_set)
  
  if(is.null(max.length)) {
    max.length <- max(nchar(input.sequences))
  }
  
  print("Padding sequences...")
  padded_sequences <- .padded.strings(input.seuqences, 
                                      max.length,
                                      padded.token = ".")
  
  
  #TODO One Hot Encoder function
  #TODO Keras Array reshape to matrix
  
  print("Tokenizing sequences...")
  sequences_tokenized <- lapply(sequences_updated, function(seq) {
    as.vector(char_to_int[strsplit(seq, "")[[1]]])
  })
  if(is.null(max.length)) {
    max.length <- max(nchar(sequences_updated))
  }
  
  

  if(convert.to.matrix) {
    print("Preparing a tokenized matrix...")
    sequences_matrix <- t(sapply(sequences_tokenized, function(x) x))
    return(sequences_matrix)
  } else {
    return(one.hot_sequences)
  }
}




