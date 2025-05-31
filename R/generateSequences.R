#' Randomly Generate Amino Acid Sequences
#' 
#' Use this to make synthetic amino acid sequences
#' for purposes of testing code, training models,
#' or noise.
#' 
#' @examples
#' generateSequences(prefix.motif = "CAS",
#'                   suffix.motif = "YF",
#'                   number.of.sequences = 100,
#'                   min.length = 8,
#'                   max.length = 16)
#'                         
#' @param prefix.motif Add defined amino acid/nucleotide sequence to the start of the generated sequences.
#' @param suffix.motif Add defined amino acid/nucleotide sequence to the end of the generated sequences
#' @param number.of.sequences Number of sequences to generate
#' @param min.length Minimum length of the final sequence. The min.length may be adjusted if 
#' incongruent with prefix.motif/suffix.motif lengths
#' @param max.length Maximum length of the final sequence
#' @param sequence.dictionary The letters to use in sequence generation (default are all amino acids)
#' 
#' @export
#' 
#' @return A vector of generated sequences
generateSequences <- function(prefix.motif = NULL,
                              suffix.motif = NULL,
                              number.of.sequences = 100,
                              min.length = 1,
                              max.length = 10,
                              sequence.dictionary = amino.acids) {
  
  # Preflight checks-----------------------------------------------------------
  prefix.motif <- if (is.null(prefix.motif)) "" else prefix.motif
  suffix.motif <- if (is.null(suffix.motif)) "" else suffix.motif
  stopifnot(is.character(prefix.motif), is.character(suffix.motif))
  stopifnot(length(number.of.sequences) == 1L, number.of.sequences >= 0L)
  
  # input sanitization
  dict <- unique(as.character(sequence.dictionary))
  bad  <- sprintf("[^%s]", paste(dict, collapse = ""))
  if (grepl(bad, prefix.motif) || grepl(bad, suffix.motif))
    stop("Prefix/suffix contain letters not in `sequence.dictionary`.")
  
  # handline length of sequences
  motif.len <- nchar(prefix.motif) + nchar(suffix.motif)
  if (motif.len > max.length)
    stop("Motifs longer than `max.length` – adjust your arguments.")
  
  min.length <- max(min.length, motif.len)
  range.len  <- min.length:max.length
  if (!length(range.len))
    stop("`min.length` exceeds `max.length` after motif adjustment.")
  
  # Generate the random sequences
  lens  <- sample(range.len, number.of.sequences, replace = TRUE)
  total <- sum(lens - motif.len)                        # pure random part
  
  rand_chars <- sample(dict, total, replace = TRUE)
  split_idx  <- rep(seq_along(lens), lens - motif.len)
  pieces     <- split.default(rand_chars, split_idx)
  
  sequences <- sprintf("%s%s%s",
                       prefix.motif,
                       vapply(pieces, paste0, collapse = "", FUN.VALUE = ""),
                       suffix.motif)
  
  sequences
}
