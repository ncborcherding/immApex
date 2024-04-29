#' Randomly Mutate Sequences of Amino Acids
#' 
#' Use this to mutate amino acid sequences
#' for purposes of testing code, training models,
#' or noise.
#' 
#' @examples
#' sequences <- generate.sequences(prefix.motif = "CAS",
#'                                 suffix.motif = "YF",
#'                                 number.of.sequences = 100,
#'                                 min.length = 8,
#'                                 max.length = 16)
#' mutated_sequences <- mutate.sequence(sequences, 
#'                                      n.sequence = 1,
#'                                      position.start = 3,
#'                                      position.end = 8)
#' 
#' 
#' @param input.sequences The amino acid sequences to use
#' @param n.sequences The number of mutated sequences to return
#' @param mutation.rate The rate of mutations to introduce into sequences
#' @param position.start The starting position to mutate along the sequence. 
#' Default NULL will start the random mutations at position 1.
#' @param position.end The ending position to mutate along the sequence.
#' Default NULL will end the random mutations at the last position.
#'@param sequence.dictionary The letters to use in sequence mutation
#' (default are all amino acids)
#' 
#' @export
#' 
#' @return A vector of mutated sequences


mutate.sequences <- function(input.sequences, 
                             n.sequences = 1, 
                             mutation.rate = 0.01,
                             position.start = NULL,
                             position.end = NULL,
                             sequence.dictionary = amino.acids[1:20]) {
  
  lapply(sequences, function(x) {
    lapply(seq_len(n.sequences), function(y) {
      .mutate.sequence(x, 
                      mutation.rate = mutation.rate,
                      position.start = position.start,
                      position.end = position.end,
                      sequence.dictionary = sequence.dictionary)
    }) -> n.mutated.sequences
    n.mutated.sequences <- unlist(n.mutated.sequences)
  }) -> output.sequences
  output.sequences <- unlist(output.sequences)
  
  return(output.sequences)
}
  

.mutate.sequence <- function(sequence, 
                             mutation.rate = mutation.rate,
                             position.start = position.start,
                             position.end = position.end,
                             sequence.dictionary = sequence.dictionary) {
    amino_acids <- strsplit(sequence, "")[[1]]
    num_mutations <- ceiling(length(amino_acids) * mutation.rate)
    if (is.null(position.start | position.start %!in% seq_len(length(amino_acids)))) {
        position.start <- 1
    }
    if (is.null(position.end) | position.end %!in% seq_len(length(amino_acids))) {
      position.end <- length(amino_acids)
    }
    
    positions_to_mutate <- sample(position.start:position.end, num_mutations)
    
    for (pos in positions_to_mutate) {
      # Select a random amino acid that is different from the current one
      possible_mutations <- setdiff(sequence.dictionary, amino_acids[pos]) # Assuming Amino acids are represented by first 24 letters
      amino_acids[pos] <- sample(possible_mutations, 1)
    }
    
    return(paste0(amino_acids, collapse = ""))
  }