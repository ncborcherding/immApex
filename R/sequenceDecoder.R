#' One Hot Decoder from One Hot Encoded Matrix or 3D Array
#' 
#' Use this to transform one hot encoded sequences back
#' into amino acid or nucleotide sequences.
#' 
#' @examples
#' decoded.sequences <- sequenceDecoder(sequence.matrix,
#'                                      sequence.dictionary = amino.acids[1:20],
#'                                      padding.symbol = ".")
#' 
#' @param sequence.matrix The encoded sequences to decode in an array opr matrix
#' @param encoder The method to prepare the sequencing information - 
#' "onehotEncoder" or "propertyEncoder"
#' @param aa.method.to.use The method or approach to use for the conversion:
#' \itemize{
#'   \item{Individual sets: atchleyFactors, crucianiProperties, FASGAI, kideraFactors, MSWHIM,
#'   ProtFP, stScales, tScales, VHSE, zScales"}
#'   \item{Multiple Sets: c("atchleyFactors", "VHSE") }
#' } 
#' @param sequence.dictionary The letters to use in sequence generation 
#' (default are all amino acids). 
#' @param padding.symbol Symbol to use for padding at the end of sequences
#' @param motif.length The length of the amino acid residues encoded 
#' (default is 1 for single amino acid encodings)
#' 
#' @export
#' @return Decoded amino acid or nucleotide sequences
sequenceDecoder <- function(sequence.matrix,
                            encoder = "onehotEncoder",
                            aa.method.to.use = NULL,
                            call.threshold = 0.5,
                            sequence.dictionary = amino.acids[1:20],
                            padding.symbol = ".") {
  
  if(encoder %!in% c("onehotEncoder", "propertyEncoder")) {
    stop("Invalid encoder provided, please select either 'onehotEncoder' or 'propertyEncoder'.")
  }
  if(encoder == "onehotEncoder") {
    decoded_sequences <- .onehotDecoder(sequence.matrix, 
                                        sequence.dictionary = sequence.dictionary,
                                        padding.symbol = padding.symbol)
  } else if (encoder == "propertyEncoder") {
    if(any(method.to.use %!in% names(apex_AA.data))) {
      stop(paste0("Please select one of the following for aa.method.to.use: ", paste(sort(names(apex_AA.data)), collapse = ", ")))
    }
    decoded_sequences <- propertyDecoder(input.sequences, 
                                         sequence.dictionary = aa.method.to.use,
                                         padding.symbol = padding.symbol)
  }
  
  
  
  return(decoded_sequences)
}

#TODO Add internal propertyDecoder
#TODO Add testthat 

#' @importFrom keras array_reshape
.onehotDecoder <- function(sequence.matrix, 
                           sequence.dictionary = amino.acids[1:20],
                           padding.symbol = padding.symbol, 
                           call.threshold = call.threshold) {
  if (inherits(sequence.matrix, "matrix")) {
    num_sequences <- nrow(sequence.matrix)
    sequence_length <- ncol(sequence.matrix) / (length(sequence.dictionary) + 1)
    sequence.matrix <- array_reshape(sequence.matrix, 
                                     c(num_sequences, sequence_length, (length(sequence.dictionary) + 1)))
  } else {
    num_sequences <- dim(sequence.matrix)[1]
    sequence_length <- dim(sequence.matrix)[2]
  }
  
  char_set <- c(sequence.dictionary, padding.symbol)
  
  decoded_sequences <- character(num_sequences)
  
  
  for (i in seq_len(num_sequences)) {
    sequence <- ""
    for (j in seq_len(sequence_length)) {
      index <- which(sequence.matrix[i, j, ] == max(sequence.matrix[i, j, ]))
      if (length(index) == 1 & max(sequence.matrix[i, j, ]) >= call.threshold) {
        sequence <- paste0(sequence, char_set[index])
      } else {
        sequence <- paste0(sequence, padding.symbol)
      }
    }
    decoded_sequences[i] <- sequence
  }
}

