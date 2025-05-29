#' Geometric Encoder from Amino Acid Strings
#' 
#' Projects each input sequence into a 20-dimensional numeric space by (i) 
#' looking up the rows of a substitution matrix (BLOSUM, PAM, or user-supplied),
#' (ii) averaging those vectors, and (iii) applying a block-diagonal 2D 
#' rotation (`θ`) across the 20 axes.
#' 
#' @section Built-in substitution matrices:
#' The package ships with ten canonical 20 × 20 matrices (stored as
#' `immapex_blosum_pam`) that can be referenced through the `method`
#' argument:
#' \itemize{
#'   \item{\strong{BLOSUM}: "BLOSUM45", "BLOSUM50", "BLOSUM62",
#'         "BLOSUM80", "BLOSUM100"}
#'   \item{\strong{PAM}:    "PAM30", "PAM40", "PAM70",
#'         "PAM120", "PAM250"}
#' }
#' 
#' @param input.sequences Character vector of AA strings.
#' @param method Character key into an internal or user list
#'   (default `"BLOSUM62"`), **or** a 20×20 numeric matrix.
#' @param substitution Optional 20×20 matrix *or* named list of matrices.
#'  When supplied, skips the built-in data lookup.
#' @param theta Rotation angle in radians for every 2-D block (default `pi/3`).
#' @param verbose Print messages corresponding to the processing step
#' 
#' @examples 
#' # Synthetic input
#' new.sequences <- generateSequences(prefix.motif = "CAS",
#'                                    suffix.motif = "YF",
#'                                    number.of.sequences = 100,
#'                                    min.length = 8,
#'                                    max.length = 16)
#'
#' # Encode with the default BLOSUM62
#' emb <- geometricEncoder(new.sequences)
#'
#' # Encode with a custom matrix
#' myMat <- matrix(runif(400), 20, 20,
#'                 dimnames = list(amino.acids, NULL))
#' emb2  <- geometricEncoder(new.sequences, substitution = myMat)
#' 
#' @return A numeric matrix with `length(input.sequences)` rows and 20
#' columns. Row order follows the input vector.
#' @export
geometricEncoder <- local({
    
    # Cache rotation matrices by theta ----------------------------------------
    .rot_cache <- new.env(parent = emptyenv())
    
    get_rotation <- function(theta) {
      key <- as.character(theta)
      if (exists(key, .rot_cache, inherits = FALSE)) return(.rot_cache[[key]])
      R2  <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
      R20 <- diag(20)
      for (i in seq(1, 20, 2)) R20[i:(i + 1), i:(i + 1)] <- R2
      .rot_cache[[key]] <- R20
      R20
    }
    
    # Internal helper to fetch the chosen substitution matrix -----------------
    fetch_matrix <- function(method, subst) {
      
      if (is.matrix(subst)) {
        if (!all(dim(subst) == 20)) stop("`substitution` matrix must be 20×20.")
        return(subst)
      }
      
      src <- if (is.null(subst)) immapex_blosum_pam else subst
      if (is.list(src) && !is.null(src[[method]])) {
        mat <- src[[method]]
        if (!all(dim(mat) == 20)) stop("Substitution matrix for '", method,
                                       "' is not 20×20.")
        return(mat)
      }
      stop("Cannot find a 20×20 matrix for method '", method, "'.")
    }
    
    # Main exported function --------------------------------------------------
    function(input.sequences,
             method = "BLOSUM62",
             substitution = NULL,
             theta = pi / 3,
             verbose = TRUE) {
      
      if (!is.character(input.sequences))
        stop("`input.sequences` must be a character vector.")
      if (anyNA(input.sequences))
        stop("NA strings are not allowed.")
      
      R20 <- get_rotation(theta)
      S   <- fetch_matrix(method, substitution)
      
      # index once for speed
      aa_lookup <- setNames(seq_along(amino.acids), amino.acids)
      
      encode_one <- function(seq) {
        idx <- aa_lookup[strsplit(seq, "", fixed = TRUE)[[1L]]]
        if (anyNA(idx))
          stop("Non-canonical amino-acid in sequence: ", seq)
        M   <- S[idx, , drop = FALSE]          
        avg <- colMeans(M)                    
        drop(R20 %*% avg)                     
      }
      
      if (verbose) message("Encoding ", length(input.sequences), " sequences ...")
      t(vapply(input.sequences, encode_one, FUN.VALUE = numeric(20)))
    }
  })
