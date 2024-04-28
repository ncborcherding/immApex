#' Adding position-specific information to sequences
#' 
#' Use this calculate positional encoding for recurrent
#' neural networks using sin/cos and position information.
#' 
#'                           
#' position.info <- positional.encoder(number.of.sequences = 1000, 
#'                                     latent.dims = 64)
#'                         
#' @param number.of.sequences The number of sequences to generate 
#' position information
#' @param latent.dims The number of latent dimensions.
#' 
#' @export
#' @return A matrix of values 

positional.encoder <- function(number.of.sequences, 
                               latent.dims = NULL) {
  
  if(length(number.of.sequences) != 1) {
    number.of.sequences <- length(number.of.sequences)
  }
  
  angle_rads <- outer(seq_len(number.of.sequences), seq_len(latent.dims), 
                      FUN=function(pos, i) {
                        pos / 10000^(2 * (i %% 2))
                      })
  
  sines <- sin(angle_rads[, seq(1, latent.dims, by=2)])
  cosines <- cos(angle_rads[, seq(2, latent.dims, by=2)])
  
  pos_encoding <- matrix(0, nrow = nrow(angle_rads), ncol = ncol(angle_rads))
  pos_encoding[, seq(1, latent.dims, by=2)] <- sines
  pos_encoding[, seq(2, latent.dims, by=2)] <- cosines
  
  return(pos_encoding)
}
