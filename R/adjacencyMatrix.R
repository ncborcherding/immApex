#' Adjacency Matrix From Amino Acid or Nucleotide Sequences
#' 
#' Calculate frequency of adjacency between residues
#' along a set of biological sequences.
#' 
#' @examples
#' new.sequences <- generateSequences(prefix.motif = "CAS",
#'                                    suffix.motif = "YF",
#'                                    number.of.sequences = 100,
#'                                    min.length = 8,
#'                                    max.length = 16)
#'                           
#' adj.matrix <- adjacencyMatrix(new.sequences,
#'                               normalize = TRUE)
#'                         
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param normalize Return the values as a function of total number of 
#' residues (\strong{TRUE}) or frequencies (\strong{FALSE})
#' @param sequence.dictionary The letters to use in sequence generation 
#' (default are all amino acids)
#' @param directed Logical; if FALSE (default) the matrix is symmetrised.
#' 
#' @export
#' @return Adjacency matrix based on input.sequences.
adjacencyMatrix <- function(input.sequences,
                            normalize = TRUE,
                            sequence.dictionary = amino.acids,
                            directed  = FALSE) {

  # Preflight checks-----------------------------------------------------------
  stopifnot(is.character(input.sequences),
            is.character(sequence.dictionary),
            length(sequence.dictionary) >= 2L)

  dict <- unique(sequence.dictionary)
  dict_lookup <- setNames(seq_along(dict), dict)

  if (length(input.sequences) == 0L)
    return(matrix(0, nrow = length(dict), ncol = length(dict),
                  dimnames = list(dict, dict)))

  # Flatten all sequences -----------------------------------------------------
  all_chars <- unlist(strsplit(input.sequences, "", fixed = TRUE), 
                      use.names = FALSE)
  bad <- is.na(dict_lookup[all_chars])
  if (any(bad))
    stop("Unknown letters found: ", paste(unique(all_chars[bad]), collapse = ", "))

  idx <- dict_lookup[all_chars]  # integer vector

  if (length(idx) < 2L)
    stop("All sequences have length 1 – no adjacencies to compute.")

  # build 'from' and 'to' vectors ---------------------------------------------
  from <- idx[-length(idx)]
  to   <- idx[-1L]

  # Tabulate to sparse COO ----------------------------------------------------
  M <- length(dict)
  bins <- (to - 1L) * M + from              # linear index for M×M matrix
  counts <- tabulate(bins, nbins = M * M)
  adj_matrix <- matrix(counts, nrow = M, ncol = M, byrow = TRUE,
                dimnames = list(dict, dict))

  ## Make symmetrical ---------------------------------------------------------
  if (!directed) {
    adj_matrix <- adj_matrix + t(adj_matrix)
  }

  ## Normalize ----------------------------------------------------------------
  if (normalize) {
    total <- sum(adj_matrix)
    if (total == 0) warning("Adjacency matrix is all zeros.")
    adj_matrix <- adj_matrix / total
  }

  return(adj_matrix)
}