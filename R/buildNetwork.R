#' Build Edit Distance Network Using Symmetric Deletion Lookup
#'
#' Constructs a weighted similarity network from biological sequences using a 
#' symmetric deletion lookup strategy combined with a banded edit-distance 
#' computation. The returned igraph object contains vertices representing the 
#' input sequences and edges representing pairs of sequences whose edit distance 
#' is less than or equal to the specified threshold. The edge attribute 
#' \code{weight} stores the computed edit distance.
#'
#' This function supports both a character vector of sequences and a data frame.
#' When provided a data frame, the user can specify the column containing sequences 
#' using the \code{sequence.column} parameter. Additionally, candidate pairs can be 
#' filtered by requiring matching \code{v.gene} and/or \code{j.gene} annotations 
#' (see \code{filter.v} and \code{filter.j}). If filtering is enabled, the corresponding 
#' gene annotation columns are required.
#'
#' @param input.data A character vector of AIR sequences, or a data frame containing 
#'   sequence data.
#' @param sequence.column A character string specifying the name of the column in 
#'   \code{input.data} that contains the sequences. Default is \code{"sequence"}. 
#'   This parameter is ignored when \code{input.data} is a character vector.
#' @param threshold An integer specifying the maximum allowed edit distance. Only 
#'   pairs of sequences with an edit distance less than or equal to this value 
#'   will be connected. Default is \code{2}.
#' @param filter.v Logical indicating whether to filter candidate pairs to only 
#'   those that have matching \code{v.gene} family annotations. Default is 
#'   \code{FALSE}. When \code{TRUE}, the input data frame must contain a column 
#'   with V gene annotations, either named \code{v.gene} or determined by 
#'   \code{.get.genes.updated}.
#' @param filter.j Logical indicating whether to filter candidate pairs to only 
#'   those that have matching \code{j.gene} family annotations. Default is 
#'   \code{FALSE}. When \code{TRUE}, the input data frame must contain a column 
#'   with J gene annotations, either named \code{j.gene} or determined by 
#'   \code{.get.genes.updated}.
#' @param technology The sequencing technology employed - \strong{'TenX'}, 
#'   \strong{'Adaptive'}, or \strong{'AIRR'}.
#' @param simplify.format If applicable, remove the allelic designation (\strong{TRUE}) or
#' retain all information (\strong{FALSE})
#' @param simplify.families If applicable, remove the hyphenated designation 
#' (\strong{TRUE}) or retain all information (\strong{FALSE})
#'
#' @return An igraph object representing the AIR similarity network. Vertices 
#' contain the original sequences (and gene annotations, if available), and each 
#' edge has a \code{weight} attribute corresponding to the computed edit distance. 
#' If no edges meet the threshold, an igraph object with only vertices is returned.
#'
#' @details
#' The function first calls a C++ routine (via Rcpp) to perform a symmetric deletion 
#' lookup, generating candidate pairs of sequences that might be within the specified 
#' edit distance. It then uses a banded dynamic programming algorithm (also implemented 
#' in C++) to compute the exact edit distance for each candidate pair. When using a 
#' data frame input, the candidate pairs can be further filtered by requiring that 
#' sequences have matching \code{v.gene} and/or \code{j.gene} values. Note that gene 
#' filtering is only applied if the corresponding filtering flag is set to \code{TRUE}.
#'
#' @examples
#' # Using a character vector of sequences:
#' sequences <- c("CASSLGTDTQYF", "CASSPGTDTQYF", "CASSLGNDTQYF", "CASRLGNDTQYF")
#' g <- buildNetwork(sequences, threshold = 2)
#' plot(g)
#'
#' # Using a data frame with a custom sequence column:
#' df <- data.frame(
#'   mySeqs = c("CASSLGTDTQYF", "CASSPGTDTQYF", "CASSLGNDTQYF", "CASRLGNDTQYF"),
#'   v.gene = c("TRBV20", "TRBV20", "TRBV12", "TRBV20"),
#'   j.gene = c("TRBJ2-7", "TRBJ2-7", "TRBJ2-1", "TRBJ2-7")
#' )
#' g_df <- buildNetwork(df, threshold = 2, filter.v = TRUE, filter.j = TRUE, sequence.column = "mySeqs")
#' plot(g_df)
#'
#' @importFrom igraph graph_from_data_frame make_empty_graph V E
#' @importFrom Rcpp evalCpp
#' @export
buildNetwork <- function(input.data, 
                         sequence.column = "sequence",
                         threshold = 2, 
                         filter.v = FALSE, 
                         filter.j = FALSE, 
                         technology = NULL, 
                         simplify.format = TRUE, 
                         simplify.families = TRUE) {
  # Handle input based on its type.
  if (is.data.frame(input.data)) {
    req_cols <- sequence.column
    
    # Determine the V gene header.
    if (filter.v) {
      v.genes.header <- .get.genes.updated(input.data, technology, "v")
      req_cols <- c(req_cols, v.genes.header)
    } else if ("v.gene" %in% colnames(input.data)) {
      v.genes.header <- "v.gene"
    } else {
      v.genes.header <- NULL
    }
    
    # Determine the J gene header.
    if (filter.j) {
      j.genes.header <- .get.genes.updated(input.data, technology, "j")
      req_cols <- c(req_cols, j.genes.header)
    } else if ("j.gene" %in% colnames(input.data)) {
      j.genes.header <- "j.gene"
    } else {
      j.genes.header <- NULL
    }
    
    # Check for required columns.
    if (!all(req_cols %in% colnames(input.data))) {
      stop(paste("Data frame must contain the following column(s):", 
                 paste(req_cols, collapse = ", ")))
    }
    
    sequences <- as.character(input.data[[sequence.column]])
    if (!is.null(v.genes.header)) {
      v_genes <- as.character(input.data[[v.genes.header]])
    } else {
      v_genes <- rep(NA, length(sequences))
    }
    if (!is.null(j.genes.header)) {
      j_genes <- as.character(input.data[[j.genes.header]])
    } else {
      j_genes <- rep(NA, length(sequences))
    }
    
  } else if (is.character(input.data)) {
    sequences <- input.data
    v_genes <- rep(NA, length(sequences))
    j_genes <- rep(NA, length(sequences))
  } else {
    stop("Input must be either a character vector or a data frame with the required columns.")
  }
  
  n <- length(sequences)
  
  if(simplify.format) {
    v_genes <- str_split(v_genes, "[*]", simplify = TRUE)[,1]
    j_genes <- str_split(j_genes, "[*]", simplify = TRUE)[,1]
  }
  
  if(simplify.familyt) {
    v_genes <- str_split(v_genes, "[-]", simplify = TRUE)[,1]
    j_genes <- str_split(j_genes, "[-]", simplify = TRUE)[,1]
  }
  
  # Get candidate pairs using the symmetric deletion lookup (implemented in C++).
  candidate_pairs <- symmetric_deletion_lookup_cpp(sequences, threshold)
  
  # Post-filter candidates to verify edit distances and obtain edge weights.
  edge_df <- post_filter_candidates(candidate_pairs, 
                                    sequences, 
                                    v_genes, 
                                    j_genes,
                                    threshold, 
                                    filter.v, 
                                    filter.j)
  
  # Build a vertices data frame including all sequences.
  vertices_df <- data.frame(name = as.character(seq_len(n)), 
                            sequence = sequences,
                            stringsAsFactors = FALSE)
  if (!all(is.na(v_genes))) {
    vertices_df$v.gene <- v_genes
  }
  if (!all(is.na(j_genes))) {
    vertices_df$j.gene <- j_genes
  }
  
  # Create the igraph object with weighted edges.
  if (nrow(edge_df) == 0) {
    g <- make_empty_graph(n)
    V(g)$sequence <- sequences
    if (!all(is.na(v_genes))) V(g)$v.gene <- v_genes
    if (!all(is.na(j_genes))) V(g)$j.gene <- j_genes
  } else {
    g <- graph_from_data_frame(d = as.data.frame(edge_df), directed = FALSE, vertices = vertices_df)
    E(g)$weight <- as.numeric(edge_df$weight)
  }
  
  return(g)
}
