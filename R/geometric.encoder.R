#' Geometric Encoder from Amino Acid Strings
#' 
#' Use this to transform amino acid sequences into a geometric 
#' encoding of the sequence.
#' 
#' @examples
#' new.sequences <- generate.sequences(prefix.motif = "CAS",
#'                                     suffix.motif = "YF",
#'                                     number.of.sequences = 100,
#'                                     min.length = 8,
#'                                     max.length = 16)
#'                           
#' sequence.matrix <- one.hot.encoder(new.sequences, 
#'                                    convert.to.matrix = TRUE)
#'                         
#' @param input.sequences The amino acid sequences to use
#' @param method.to.use The method or approach to use for the conversion: 
#' "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", 
#' "PAM40", "PAM70", "PAM120", or "PAM250"  
#' @param theta angle to use for geometric transformation
#' 
#' @export
#' @return Geometric encoded amino acid sequences in a matrix


geometric.encoder <- function(input.sequences, 
                               method.to.use = "BLOSUM62",
                               theta = pi) {
  possible.methods <- c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", 
                        "BLOSUM100", "PAM30","PAM40", "PAM70", "PAM120", "PAM250")
  
  if (method.to.use %!in% possible.methods) {
    stop(paste0("Please select a method.to.us from the following options: ", 
                paste(possible.methods, collapse = ", ")))
  }
  
  if(any(unlist(strsplit(input.sequences[1:10], "")) %!in% amino.acids[1:20])) {
    stop("geometric.encoder() works only on the sequences with the conventional 20 amino acids")
  }
  
  dim <- 20
  rotation_matrix <- diag(dim)
  
  # Generate ten 2D rotation matrices and insert them into the 20D rotation matrix
  rotation_2d = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol=2)
  
  for (i in seq(1, dim, by=2)) {
    rotation_matrix[i:(i+1), i:(i+1)] = rotation_2d
  }
  
  print("Performing geometric transformation...")
  geometric_sequences <- lapply(input.sequences, function(x) {
    if(is.na(x)) {
      tmp <- rep(0, 20)
    } else {
      tmp <- .encode_sequence(x, 
                              dim = dim, 
                              rotation_matrix = rotation_matrix,
                              method = method.to.use)
    }
    tmp
  })
  
  geommetric_matrix <- do.call(rbind,geometric_sequences)
  return(geommetric_matrix)
}

#Create BLOSUM/PAM Matrix
.aa_to_submatrix <- function(sequence,
                             method) {
  aa_indices <- match(strsplit(as.character(sequence), '')[[1]], amino.acids[1:20])
  return(apex_blosum.pam.matrices[[method]][aa_indices, ])
}

# Main function to encode each sequence
#' @importFrom matrixStats colWeightedMeans
.encode_sequence <- function(sequence, 
                             dim, 
                             rotation_matrix,
                             method) {
  
  blosum_vectors <- .aa_to_submatrix(sequence,
                                     method)
  
  # Initialize an empty matrix to store the transformed points
  transformed_points <- matrix(nrow=0, ncol=dim)
  
  # Apply the unitary transformation
  for (i in 1:(dim(blosum_vectors)[1])) {
    point <- t(blosum_vectors[i, ])
    
    # Apply the rotation matrix
    transformed_point <- rotation_matrix %*% matrix(point, nrow=dim, ncol=1)
    
    # Collect the transformed points
    transformed_points <- rbind(transformed_points, t(transformed_point))
  }
  
  # Average the transformed points to get a single 20D point
  avg_point <- colWeightedMeans(transformed_points)
  
  return(avg_point)
}

