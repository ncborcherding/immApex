#' Universal Amino-acid Sequence Encoder 
#'
#' `sequenceEncoder()` is a high-level R wrapper around the low-level
#' `encodeSequences_cpp()` engine (see `?encodeSequences_cpp`).  It converts a
#' character vector of amino-acid sequences into either  
#' (1) **one-hot** representations, or  
#' (2) **property-based** representations (e.g. hydropathy, charge, Atchley
#' factors)  
#'
#' The function performs basic argument validation, constructs (or accepts) a
#' property matrix when `mode = "property"`, and then calls the C++ back-end.
#' The resulting list contains a 3-D array (`cube`) and a flattened 2-D matrix
#' (`flattened`) that can feed directly into downstream machine-learning
#' pipelines (autoencoders, random forests, etc.).
#'
#' @section Property mode:
#' If you supply `property.matrix` directly, it **must** be a numeric matrix
#' whose **rows correspond to the 20 canonical amino acids in the order of
#' `sequence.dictionary`** and whose columns are the property scales.  
#'
#' @param input.sequences `character` vector. Sequences (uppercase
#' single-letter code).  Gaps or unknown symbols are replaced by `padding.symbol`.
#' @param mode Either `"onehot"` (default) or `"property"`.
#' @param property.set *Optional* `character` vector of property names to
#' extract from \pkg{Peptides} (ignored in `"onehot"` mode).  Ignored if 
#' `property.matrix` is supplied.
#' @param property.matrix *Optional* numeric matrix (`20 × P`). Overrides 
#' `property.set`.
#' @param sequence.dictionary Character vector of the nucleotide or amino-acid 
#' alphabet (default = 20 standard residues). The order defines the row order of 
#' `property.matrix`.
#' @param padding.symbol Single character used for right-padding. Must not be
#'   one of the sequence.dictionary.
#' @param summary.fun For property mode only: currently `"mean"` or `""`
#' (empty string) for no summary.  When `"mean"`, the computes the per-sequence mean across
#' positions and returns it as an extra element `summary`.
#' @param max.length  Positive integer. Sequences longer than this are
#' truncated; shorter ones are right-padded with `padding.symbol`. If `NULL` 
#' (default) the longest input sequence sets the maximum.
#' @param nthreads Number of threads to request (passed to the C++ back-end).
#' Honoured only when the package is compiled with OpenMP; otherwise silently 
#' coerced to `1`.  Default = `parallel::detectCores()`.
#' @param verbose Logical. If `TRUE` (default) prints a short progress message.
#'
#' @return A named `list` with the following elements: 
#' \describe{
#'   \item{`cube`}{Numeric array with dimensions  
#'                 `c(depth, max.length, length(sequences))`.  
#'                 Depth = *sequence.dictionary size* in `"onehot"` mode, or
#'                 *number of property scales* in `"property"` mode.}
#'   \item{`flattened`}{Numeric matrix, `length(sequences)` rows by
#'                      `depth × max.length` columns.}
#'   \item{`summary`}{(Only when `summary.fun != ""`) Numeric matrix,
#'                   `length(sequences) × depth`, containing the requested
#'                   statistic.}
#'   \item{`sequence.dictionary`, `padding.symbol`, `threads`}{Metadata echoed from the call.}
#' }
#'
#'
#' @examples
#' aa <- c("CARDRST", "YYYGMD", "ACACACAC")
#' enc <- sequenceEncoder(aa, 
#'                        mode = "onehot")
#'
#' prop <- sequenceEncoder(aa,
#'                         mode = "property",
#'                         property.set = c("Atchley"),
#'                         summary.fun  = "mean")
#'
#' @export
sequenceEncoder <- function(input.sequences,
                            mode             = c("onehot", "property"),
                            property.set     = NULL,
                            property.matrix  = NULL,
                            sequence.dictionary = amino.acids,       
                            padding.symbol   = ".",
                            summary.fun      = "",
                            max.length       = NULL,
                            nthreads         = parallel::detectCores(),
                            verbose          = TRUE) {
  
  mode <- match.arg(mode)
  
  #  Property matrix construction / validation --------------------------------
  if (mode == "property") {
    if (!is.null(property.set)) {
      prop_mat <- t(.aa.property.matrix(property.set))
    } else if (!is.null(property.matrix)) {
      if (!is.matrix(property.matrix) || !is.numeric(property.matrix) ||
          nrow(property.matrix) != length(sequence.dictionary)) {
        stop(sprintf(
          "`property.matrix` rows (%d) must match `sequence.dictionary` length (%d).",
          nrow(property.matrix), length(sequence.dictionary)
        ))
      }
      prop_mat <- property.matrix
    } else {
      stop("In `property` mode, you must supply either `property.set` 
           (as a character vector) or `property.matrix` (as a numeric matrix).")
    }
  } else {
    prop_mat <- NULL
  }
  
  # Max length set -----------------------------------------------------------
  if (is.null(max.length))
    max.length <- max(nchar(input.sequences), 1L)
  
  # Encoding sequences --------------------------------------------------------
  if (verbose)
    message(sprintf("[sequenceEncoder] Encoding %d sequence%s (%s mode, L≤%d)…",
                    length(input.sequences), ifelse(length(input.sequences) == 1L, "", "s"),
                    mode, max.length))
  
  out <- encodeSequences_cpp(
    sequences      = input.sequences,
    mode           = mode,
    alphabet       = sequence.dictionary,
    prop_mat_      = prop_mat,
    pad_token      = padding.symbol,
    summary        = summary.fun,
    max_len        = max.length,
    nthreads       = nthreads
  )
  
  return(out)
}

#' @rdname sequenceEncoder                                   
#' @aliases onehotEncoder                                     
#' @export
onehotEncoder <- function(...,
                          mode = "onehot") {    
  sequenceEncoder(..., mode = "onehot")
}

#' @rdname sequenceEncoder
#' @aliases propertyEncoder
#' @export
propertyEncoder <- function(...,
                            mode = "property") {
  sequenceEncoder(..., mode = "property")
}