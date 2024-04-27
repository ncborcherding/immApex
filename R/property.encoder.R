"crucianiProperties" = crucianiProperties,
"fasgaiVectors" = fasgaiVectors,
"kideraFactors" = kideraFactors,
"mswhimScores" = mswhimScores,
"protFP" protFP,
"stScales" = stScales,
"tScales" = tScales,
"vhseScales" = vhseScales, 
"zScales" = zScales

property.encoder <- function(input.sequences, 
                             max.length = NULL,
                             method.to.use = NULL,
                             convert.to.matrix = TRUE) {
  
  propertyFunc <- switch(method.to.use, 
                         "crucianiProperties" = crucianiProperties,
                         "fasgaiVectors" = fasgaiVectors,
                         "kideraFactors" = kideraFactors,
                         "mswhimScores" = mswhimScores,
                         "protFP" protFP,
                         "stScales" = stScales,
                         "tScales" = tScales,
                         "vhseScales" = vhseScales, 
                         "zScales" = zScales)
  
  #TODO Add Apex AA list
  #TODO Extract Vectors
  #TODO Normalize Vectors
  vectors <- AAdata$FASGAI
  
  if(is.null(max.length)) {
    max.length <- max(nchar(input.sequences))
  }
  
  #How to pad motifs
  print("Padding sequences...")
  padded_sequences <- .padded.strings(strings = input.sequences, 
                                      max.length = max.length,
                                      padded.token = ".",
                                      concatenate = TRUE)
  
  print("Property-based Encoding sequences...")
  property_sequences <- .convert.property(unlist(padded_sequences),
                                          max.length = max.length,
                                          vectors = vectors)
  
  if(convert.to.matrix) {
    print("Preparing a matrix...")
    property_matrix <- array_reshape(property_sequences, c(dim(property_sequences)[1], dim(property_sequences)[2]*dim(property_sequences)[3]))
    return(property_matrix)
  } else {
    return(property_sequences)
  }
  
}

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
  
  dimnames(property_array) <- list(paste0("Seq_", 1:length(sequences)),
                                  paste0("Pos_", 1:max.length),
                                  names(vectors))
  return(property_array)
}
