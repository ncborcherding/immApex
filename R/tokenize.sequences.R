#' Generate Tokenized Sequences from Amino Acid String
#' 
#' Use this to transform amino acid sequences into 
#' tokens in preparing for deep learning models. 
#' 
#' @examples
#' new.sequences <- generate.sequences(prefix.motif = "CAS",
#'                                     suffix.motif = "YF",
#'                                     number.of.sequences = 100,
#'                                     min.length = 8,
#'                                     max.length = 16)
#'                           
#'sequence.matrix <- tokenize.sequences(new.sequences, 
#'                                      add.startstop = TRUE,
#'                                      start.token = "!",
#'                                      stop.token = "^", 
#'                                      convert.to.matrix = TRUE)
#'                         
#' @param input.sequences The amino acid sequences to use
#' @param add.startstop Add start and stop tokens to the sequence
#' @param start.token The character to use for the start token
#' @param stop.token The character to use for the stop token
#' @param max.length Additional length to pad, NULL will pad sequences 
#' to the max length of input.sequences
#' @param convert.to.matrix Return a matrix (TRUE) or a vector (FALSE)
#' 
#' @export
#' @return Tokenized amino acid sequences in a matrix or vector
tokenize.sequences <- function(input.sequences, 
                               add.startstop = TRUE,
                               start.token = "!",
                               stop.token = "^", 
                               max.length = NULL,
                               convert.to.matrix = TRUE) {
  
  if(add.startstop) {
    char_set <- c(start.token,amino.acids, stop.token)
    print("Adding start and stop tokens...")
    sequences_updated <- sapply(input.sequences, function(x) paste(start.token, x, stop.token, sep = ""))
  } else {
    char_set <- c(amino.acids)
  }
  
  # Create a mapping of amino acids to integers
  char_to_int <- setNames(seq_along(char_set), char_set)
  
  print("Tokenizing sequences...")
  sequences_tokenized <- lapply(sequences_updated, function(seq) {
    as.vector(char_to_int[strsplit(seq, "")[[1]]])
  })
  if(is.null(max.length)) {
    max.length <- max(nchar(sequences_updated))
  }
  
  print("Padding sequences...")
  sequences_tokenized <- .padded.strings(sequences_tokenized, 
                                         max.length,
                                         padded.token = length(char_to_int) + 1,
                                         concatenate = FALSE)
  
  if(convert.to.matrix) {
    print("Preparing a tokenized matrix...")
    sequences_matrix <- t(sapply(sequences_tokenized, function(x) x))
    return(sequences_matrix)
  } else {
    return(sequences_tokenized)
  }
}




