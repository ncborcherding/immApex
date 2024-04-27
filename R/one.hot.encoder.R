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
#' @importFrom keras array_reshape
#' 
#' @export
#' @return One hot encoded amino acid sequences in a matrix or 3D array
one.hot.encoder <- function(input.sequences, 
                            max.length = NULL,
                            split.length = 1,
                            convert.to.matrix = TRUE,
                            sequence.dictionary = amino.acids[1:20]) {
  if(split.length == 1) {
    char_set <- c(sequence.dictionary, ".")
  } else {
    all_motifs <- expand.grid(replicate(split.length, sequence.dictionary, simplify = FALSE))
    unique_motifs <- unique(apply(all_motifs, 1, paste, collapse = ""))
    char_set <- c(unique_motifs, ".")
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
                                      padded.token = ".",
                                      concatenate = TRUE)
  
  #How to convert to motifs
  print("One Hot Encoding sequences...")
  onehot_sequences <- .convert.one.hot(unlist(padded_sequences),
                                       max.length = max.length,
                                       char_set = char_set)

  if(convert.to.matrix) {
    print("Preparing a matrix...")
    onehot_matrix <- array_reshape(onehot_sequences, c(dim(onehot_sequences)[1], dim(onehot_sequences)[2]*dim(onehot_sequences)[3]))
    return(onehot_matrix)
  } else {
    return(onehot_sequences)
  }
}

#TODO Allow for motif or single AA encoding



.convert.one.hot <- function(sequences, 
                             split.length = 1,
                             max.length,
                             char_set = NULL) {
  
  one_hot_array <- array(0, dim = c(length(sequences), max.length, length(char_set)))
  for (i in seq_len(length(sequences))) {
    chars <- strsplit(sequences[i], "")[[1]]
    valid_indices <- match(chars, char_set, nomatch = length(char_set) + 1) #NoMatch will not be recorded
    for(t in seq_along(chars)) {
      one_hot_array[i, t, valid_indices[t]] <- 1
    }
  }
  
  dimnames(one_hot_array) <- list(paste0("Seq_", 1:length(sequences)),
                                  paste0("Pos_", 1:max.length),
                                  c(char_set))
  return(one_hot_array)
}
