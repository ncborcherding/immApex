#' Geometric Encoder from Amino Acid Strings
#' 
#' #' Projects sequences into a 20D space using a substitution matrix,
#' vector averaging, and rotation. This version is optimized for speed
#' by using fully vectorized operations.
#'
#' @param input.sequences Character vector of AA strings.
#' @param method Character key for a built-in substitution matrix 
#'   (e.g., "BLOSUM62"), or a 20x20 numeric matrix itself.
#' @param substitution Optional named list of matrices to look up `method` in.
#'   If `method` is a matrix, this is ignored.
#' @param theta Rotation angle in radians (default `pi/3`).
#'
#' @return A numeric matrix with `length(input.sequences)` rows and 20 columns.
#' @export
geometricEncoder <- function(input.sequences,
                             method = "BLOSUM62",
                             substitution = NULL,
                             theta = pi / 3) {
  
  # Preflight checks-----------------------------------------------------------
  if (!is.character(input.sequences)) {
    stop("`input.sequences` must be a character vector.")
  }
  if (length(input.sequences) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 20))
  }
  if (anyNA(input.sequences) || any(nchar(input.sequences) == 0)) {
    stop("NA or empty strings are not allowed in `input.sequences`.")
  }
  
  # 2. Fetch Substitution Matrix 
  fetch_matrix <- function(m, s) {
    if (is.matrix(m)) {
      if (!all(dim(m) == 20) || is.null(rownames(m))) {
        stop("If `method` is a matrix, it must be 20x20 with amino acid rownames.")
      }
      return(m)
    }
    if (!exists("immapex_blosum.pam.matrices")) {
      stop("Built-in matrix data `immapex_blosum.pam.matrices` not found.")
    }
    
    src <- if (is.null(s)) immapex_blosum.pam.matrices else s
    
    mat <- src[[m]]
    if (is.null(mat)) stop("Cannot find matrix for method '", m, "'.")
    if (!all(dim(mat) == 20) || is.null(rownames(mat))) {
      stop("Matrix for '", m, "' is not a 20x20 matrix with rownames.")
    }
    mat
  }
  
  S <- fetch_matrix(method, substitution)
  aa_lookup <- setNames(seq_len(nrow(S)), rownames(S))
  
  # 3. Vectorized Averaging 
  # Get lengths of each sequence
  seq_lengths <- nchar(input.sequences)
  
  # Create a grouping factor to identify which sequence each amino acid belongs to
  group_id <- rep.int(seq_along(seq_lengths), seq_lengths)
  
  # Ungroup all sequences into a single character vector
  all_chars <- unlist(strsplit(input.sequences, "", fixed = TRUE), use.names = FALSE)
  
  # Perform lookup for all amino acids at once
  all_indices <- aa_lookup[all_chars]
  
  # Check for any non-canonical amino acids
  if (anyNA(all_indices)) {
    bad_chars <- unique(all_chars[is.na(all_indices)])
    stop("Non-canonical amino acid(s) found: ", paste(bad_chars, collapse = ", "))
  }
  
  # Get all corresponding vectors from the substitution matrix
  all_vectors <- S[all_indices, , drop = FALSE]
  summed_vectors <- rowsum(all_vectors, group_id, reorder = FALSE)
  avg_vectors <- summed_vectors / seq_lengths
  
  # 4. Vectorized Rotation 
  R2 <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  R20 <- diag(20)
  for (i in seq(1, 19, 2)) {
    R20[i:(i + 1), i:(i + 1)] <- R2
  }
  rotated_vectors <- avg_vectors %*% t(R20)
  
  return(rotated_vectors)
}