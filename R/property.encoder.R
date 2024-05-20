#' Encoder from Amino Acid String by Properties
#' 
#' Use this to transform amino acid sequences a
#' a matrix by amino acid properties derived from 
#' dimensional reduction strategies
#' 
#' @examples
#' new.sequences <- generate.sequences(prefix.motif = "CAS",
#'                                     suffix.motif = "YF",
#'                                     number.of.sequences = 100,
#'                                     min.length = 8,
#'                                     max.length = 16)
#'                           
#' sequence.matrix <- property.encoder(new.sequences, 
#'                                     method.to.use = "VHSE",
#'                                     convert.to.matrix = TRUE)
#'                         
#' @param input.sequences The amino acid sequences to use
#' @param max.length Additional length to pad, NULL will pad sequences 
#' to the max length of input.sequences
#' @param method.to.use The method or approach to use for the conversion: 
#' "atchleyFactors", "crucianiProperties", "FASGAI", "kideraFactors", "MSWHIM",
#' "ProtFP", "stScales", "tScales", "VHSE", or "zScales". Multiple embeddings are
#' possible, provide a vector of approaches - c("atchleyFactors", "VHSE")           
#' @param convert.to.matrix Return a matrix (TRUE) or a 3D array (FALSE)
#' @param summary.function Return a matrix that summarize the amino acid method/property
#' Available summaries include: "median", "mean", "sum", variance ("vars"), or 
#' Median Absolute Deviation ("mads")
#' @param padding.symbol Symbol to use for padding at the end of sequences
#' @importFrom keras array_reshape
#' 
#' @export
#' @return Converted amino acid sequences by property in a matrix or 3D array

property.encoder <- function(input.sequences, 
                             max.length = NULL,
                             method.to.use = NULL,
                             convert.to.matrix = TRUE,
                             summary.function = NULL,
                             padding.symbol = ".") {
  if(any(method.to.use %!in% names(apex_AA.data))) {
    stop(paste0("Please select one of the following for method.to.use: ", paste(sort(names(apex_AA.data)), collapse = ", ")))
  }
  vectors <- apex_AA.data[method.to.use]
  vector.names <- as.vector(unlist(lapply(vectors, names)))
  vectors <- do.call(c, vectors)
  names(vectors) <- vector.names
  
  #TODO Think about other normalization
  vectors <- lapply(vectors, .min.max.normalize)
  
  if(is.null(max.length)) {
    max.length <- max(nchar(input.sequences))
  }
  
  #How to pad motifs
  print("Padding sequences...")
  padded_sequences <- .padded.strings(strings = input.sequences, 
                                      max.length = max.length,
                                      padded.token = padding.symbol,
                                      concatenate = TRUE)
  
  print("Property-based Encoding sequences...")
  property_sequences <- .convert.property(unlist(padded_sequences),
                                          max.length = max.length,
                                          vectors = vectors)
  
  if(!is.null(summary.function)) {
    if (tolower(summary.function) %!in% c("median", "mean", "sum", "vars", "mads")) {
      stop("Please select one of the following summary.function options: 'median', 'mean', 'sum', 'vars', or 'mads'")
    }
    print(paste0("Summarising properties using ", summary.function, "..."))
    stat.function <- .get.stat.function(summary.function)
    stat_matrix <- matrix(nrow = dim(property_sequences)[1], ncol = dim(property_sequences)[3])
    
    for (i in seq_len(length(input.sequences))) {
      stat_matrix[i, ] <- stat.function(property_sequences[i, , ], na.rm = TRUE)
    }
    colnames(stat_matrix) <- dimnames(property_sequences)[[3]]
    return(stat_matrix)
  }
  
  if(convert.to.matrix) {
    print("Preparing a matrix...")
    property_matrix <- array_reshape(property_sequences, c(dim(property_sequences)[1], dim(property_sequences)[2]*dim(property_sequences)[3]))
    colnames(property_matrix) <- array.dimnamer(property_sequences)
    return(property_matrix)
  } else {
    return(property_sequences)
  }
  
}

#Converting the sequence into numerical matrix
.convert.property <- function(sequences, 
                              max.length,
                              vectors = vectors) {
  property_array <- array(0, dim = c(length(sequences), max.length, length(vectors)))
  
  for (i in seq_len(length(sequences))) {
    transformed <-sapply(names(vectors), function(scale) {
      vectors[[scale]][strsplit(sequences[[i]], "")[[1]]]
    })
    transformed[is.na(transformed)] <- 0
    property_array[i,,] <- transformed
  }
  
  dimnames(property_array) <- list(paste0("Seq.", 1:length(sequences)),
                                  paste0("Pos.", 1:max.length),
                                  names(vectors))
  return(property_array)
}
