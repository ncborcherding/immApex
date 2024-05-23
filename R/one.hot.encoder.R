#' One Hot Encoder from Amino Acid or Nucleotide Strings
#' 
#' Use this to transform amino acid or nucleotide sequences 
#' into a one hot encoding of the sequence.
#' 
#' @examples
#' new.sequences <- generate.sequences(prefix.motif = "CAS",
#'                                     suffix.motif = "YF",
#'                                     number.of.sequences = 100,
#'                                     min.length = 8,
#'                                     max.length = 16)
#'                           
#' sequence.matrix <- one.hot.encoder(new.sequences, 
#'                                    convert.to.matrix = TRUE)
#'                         
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param max.length Additional length to pad, NULL will pad sequences 
#' to the max length of input.sequences
#' @param motif.length The length of the amino acid residues to encode -
#' a motif.length = 1 produces single amino acid encodings
#' @param convert.to.matrix Return a matrix (\strong{TRUE}) or a 3D array (\strong{FALSE})
#' @param sequence.dictionary The letters to use in sequence generation 
#' (default are all amino acids). This will be overrode if using a
#'  motif approach (\strong{motif.length} > 1)
#' @param padding.symbol Symbol to use for padding at the end of sequences
#' 
#' @importFrom keras array_reshape
#' @importFrom stats setNames
#' 
#' @export
#' @return One hot encoded sequences in a matrix or 3D array
one.hot.encoder <- function(input.sequences, 
                            max.length = NULL,
                            motif.length = 1,
                            convert.to.matrix = TRUE,
                            sequence.dictionary = amino.acids[1:20],
                            padding.symbol = ".") {
  
  char_set <- c(sequence.dictionary, padding.symbol)
  
  if(.check.sequences(input.sequences, sequence.dictionary)) {
    stop("The sequence.dictionary does not cover the input.sequences, please modify or add the additional characters.")
  }

  if (motif.length > 1) {
    all_motifs <- expand.grid(replicate(motif.length, char_set, simplify = FALSE))
    char_set <- unique(apply(all_motifs, 1, paste, collapse = ""))
  }
  # Create a mapping of amino acids to integers
  char_to_int <- setNames(seq_along(char_set), char_set)
  
  if(is.null(max.length)) {
    max.length <- max(nchar(input.sequences))
  }
  
  #How to pad motifs
  print("Padding sequences...")
  padded_sequences <- .padded.strings(strings = input.sequences, 
                                      max.length = max.length,
                                      padded.token = padding.symbol,
                                      concatenate = TRUE)
  
  #How to convert to motifs
  print("One Hot Encoding sequences...")
  onehot_sequences <- .convert.one.hot(unlist(padded_sequences),
                                       max.length = max.length,
                                       motif.length = motif.length,
                                       char_set = char_set)

  if(convert.to.matrix) {
    print("Preparing a matrix...")
    onehot_matrix <- array_reshape(onehot_sequences, c(dim(onehot_sequences)[1], dim(onehot_sequences)[2]*dim(onehot_sequences)[3]))
    colnames(onehot_matrix) <- array.dimnamer(onehot_sequences)
    return(onehot_matrix)
  } else {
    return(onehot_sequences)
  }
}



.convert.one.hot <- function(sequences, motif.length = 1, max.length, char_set = NULL) {
  # Initialize the one-hot array with appropriate dimensions
  one_hot_array <- array(0, dim = c(length(sequences), max.length, length(char_set)))
  
  # Extract all subsequences from each sequence
  subsequences <- substring.extractor(sequences, motif.length)
  
  # Apply one-hot encoding
  for (i in seq_along(subsequences)) {
    chars <- subsequences[[i]]
    valid_indices <- match(chars, char_set)
    for(t in seq_along(chars)) {
      if (!is.na(valid_indices[t])) {
        one_hot_array[i, t, valid_indices[t]] <- 1
      }
    }
  }
  
  dimnames(one_hot_array) <- list(paste0("Seq.", 1:length(sequences)),
                                  paste0("Pos.", 1:max.length),
                                  char_set)
  return(one_hot_array)
}
