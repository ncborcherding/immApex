#' Relative Residue Frequencies at Every Position
#'
#' Quickly computes the per-position relative frequency of each symbol
#' (amino-acid or nucleotide) in a set of biological sequences.  Variable-length
#' strings are padded to a common width so the calculation is entirely
#' vectorised (one logical comparison + one `colSums()` per residue).
#'
#' @param sequences  Character vector of sequences (amino acid or 
#' nucleotide)
#' @param max.length Integer.  Pad/trim to this length. Defaults to
#'  `max(nchar(sequences))`.
#' @param sequence.dictionary Vector of valid residue symbols that should be
#' tracked (defaults to the 20 canonical amino acids; supply
#' `c("A","C","G","T","N")` etc. for nucleotides).
#' @param padding.symbol Single character used for right-padding. **Must not**
#' be present in `sequence.dictionary`.
#' @param tidy Logical; if `TRUE` a long-format `data.frame` is returned
#' instead of a matrix (useful for plotting with *ggplot2*).
#'
#' @return Either
#' \itemize{
#'   \item A numeric matrix of dimension `length(sequence.dictionary)` ×
#'         `max.length`, whose columns sum to 1, **or**
#'   \item A `data.frame` with columns *position*, *residue*, *frequency* when
#'         `tidy = TRUE`.
#' }
#' @examples
#' set.seed(1)
#' seqs <- c("CASSLGQGAETQYF", "CASSPGQGDYEQYF", "CASSQETQYF")
#' rel.freq <- calculateFrequency(seqs)
#' head(rel.freq[, 1:5])                  
#'
#' ## Nucleotide example
#' dna <- c("ATGCC", "ATGAC", "ATGGC")
#' calculateFrequency(dna,
#'                    sequence.dictionary = c("A","C","G","T"),
#'                    padding.symbol = "-",
#'                    tidy = TRUE)
#'
#' @export
calculateFrequency <- function(sequences,
                               max.length = NULL,
                               sequence.dictionary = amino.acids,
                               padding.symbol = ".",
                               tidy = FALSE) {
  # Preflight checks-----------------------------------------------------------
  stopifnot(is.character(sequences),
            nchar(padding.symbol) == 1L,
            padding.symbol %!in% sequence.dictionary)
  
  # 1. Pad to a rectangular character matrix 
  if (is.null(max.length))
    max.length <- max(nchar(sequences), 1L)
  
  padded <- .padded.strings(sequences,
                            max.length = max.length,
                            pad = padding.symbol,
                            collapse  = TRUE)
  
  seq_mat <- do.call(rbind, strsplit(unlist(padded), ""))  # nSeq × max.length
  nSeq <- length(padded)
  
  # 2.  Fast frequency calculation (one pass per residue) 
  res_mat <- matrix(0,
                    nrow = length(sequence.dictionary) + 1,
                    ncol = max.length,
                    dimnames = list(c(sequence.dictionary, padding.symbol),
                                    paste0("Pos.", seq_len(max.length))))
  
  for (residue in c(sequence.dictionary, padding.symbol)) {
    # logical comparison is vectorised; colSums is C-level
    res_mat[residue, ] <- colSums(seq_mat == residue) / nSeq
  }
  
  # 3.  Optional tidy reshaping 
  if (tidy) {
    res_mat <- as.data.frame(as.table(res_mat),
                             stringsAsFactors = FALSE,
                             responseName = "frequency")
    names(res_mat) <- c("residue", "position", "frequency")
    res_mat$position <- as.integer(sub("Pos\\.", "", res_mat$position))
  }
  
  res_mat
}
