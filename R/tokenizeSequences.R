#' Generate Tokenized Sequences from Amino Acid String
#' 
#' Use this to transform amino acid sequences into 
#' tokens in preparing for deep learning models. 
#' 
#' @examples
#' new.sequences <- generateSequences(prefix.motif = "CAS",
#'                                    suffix.motif = "YF",
#'                                    number.of.sequences = 100,
#'                                    min.length = 8,
#'                                    max.length = 16)
#'                           
#'sequence.matrix <- tokenizeSequences(new.sequences, 
#'                                     add.startstop = TRUE,
#'                                     start.token = "!",
#'                                     stop.token = "^", 
#'                                     convert.to.matrix = TRUE)
#'                         
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param add.startstop Add start and stop tokens to the sequence
#' @param start.token The character to use for the start token
#' @param stop.token The character to use for the stop token
#' @param max.length Additional length to pad, NULL will pad sequences 
#' to the max length of input.sequences
#' @param convert.to.matrix Return a matrix (TRUE) or a vector (FALSE)
#' @param verbose Print messages corresponding to the processing step
#' 
#' @export
#' @return Tokenize sequences in a matrix or vector
tokenizeSequences <- function(input.sequences, 
                              add.startstop = TRUE,
                              start.token = "!",
                              stop.token = "^", 
                              max.length = NULL,
                              convert.to.matrix = TRUE,
                              verbose = TRUE) {
  
  if(add.startstop) {
    char_set <- c(start.token,amino.acids, stop.token)
    message("Adding start and stop tokens...")
    sequences_updated <- vapply(input.sequences, 
                                function(x) paste(start.token, x, stop.token, sep = ""), 
                                FUN.VALUE = character(1))
  } else {
    char_set <- c(amino.acids)
    sequences_updated <- input.sequences
  }
  
  # Create a mapping of amino acids to integers
  char_to_int <- setNames(seq_along(char_set), char_set)
  
  if(verbose) {
    message("Tokenizing sequences...")
  }
  sequences_tokenized <- lapply(sequences_updated, function(seq) {
    as.vector(char_to_int[strsplit(seq, "")[[1]]])
  })
  if(is.null(max.length)) {
    max.length <- max(nchar(sequences_updated))
  }
  
  if(verbose){
    message("Padding sequences...")
  }
  sequences_tokenized <- .padded.strings(sequences_tokenized, 
                                         max.length,
                                         padded.token = length(char_to_int) + 1,
                                         concatenate = FALSE)
  
  if(convert.to.matrix) {
    if(verbose) {
      message("Preparing a tokenized matrix...")
    }
    sequences_matrix <- do.call(rbind, sequences_tokenized)
    return(sequences_matrix)
  } else {
    return(sequences_tokenized)
  }
}




