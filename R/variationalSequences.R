#' Generate Similar Sequences using Variational Autoencoder
#' 
#' Use this to simulate sequences using a variational 
#' autoencoder (VAE) and perturbation of the probability distributions.
#' 
#' @examples
#' sequences <- generateSequences(prefix.motif = "CAS",
#'                                suffix.motif = "YF",
#'                                number.of.sequences = 100,
#'                                min.length = 8,
#'                                max.length = 16)
#'                                 
#' variate_sequences <- variationalSequences(sequences, 
#'                                           encoder = "onehotEncoder",
#'                                           layers = 2, 
#'                                           hidden.dims = c(256, 128)
#'                                           latent.dim = 16,
#'                                           batch.size = 16)
#' 
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param encoder The method to prepare the sequencing information - 
#' "onehotEncoder" or "propertyEncoder"
#' @param aa.method.to.use The method or approach to use for the conversion:
#' \itemize{
#'   \item{Individual sets: atchleyFactors, crucianiProperties, FASGAI, kideraFactors, MSWHIM,
#'   ProtFP, stScales, tScales, VHSE, zScales"}
#'   \item{Multiple Sets: c("atchleyFactors", "VHSE") }
#' } 
#' @param layers The number of hidden layers to employ within the VAE
#' @param number.of.sequences Number of sequences to generate
#' @param hidden.dims A vector of the neurons to use in the hidden layers, 
#' The length needs to match the number of layers
#' @param latent.dim The size of the latent dimensions
#' @param batch.size The batch size to use for VAE training
#' @param epochs The number of epochs to use in VAE training
#' @param learning.rate The learning rate to use in VAE training
#' @param epsilon.std The epsilon to use in VAE training
#' @param null.threshold The null threshold to use in VAE training
#' @param call.threshold The relative strictness of sequence calling
#'  with higher values being more stringent
#' @param activation.function The activation for the dense connected layers
#' @param optimizer The optimizer to use in VAE training
#' @param disable.eager.execution Disable the eager execution parameter for
#' tensorflow.
#' @param sequence.dictionary The letters to use in sequence mutation
#' (default are all amino acids)
#' 
#' @importFrom keras layer_dense layer_concatenate layer_lambda 
#' keras_model compile fit k_sum layer_input loss_binary_crossentropy 
#' optimizer_adadelta optimizer_adagrad optimizer_adam optimizer_adamax 
#' optimizer_ftrl optimizer_nadam optimizer_rmsprop optimizer_sgd 
#' backend
#' @importFrom dplyr %>%
#' @importFrom stats predict
#' @importFrom tensorflow tf
#' @export 
#' @return A vector of mutated sequences



variationalSequences <-function(input.sequences,
                                number.of.sequences = 100,
                                encoder = "onehotEncoder",
                                aa.method.to.use = NULL,
                                layers = 2, 
                                hidden.dims = c(256, 128),
                                latent.dim = 16,
                                batch.size = 16,
                                epochs = 30,
                                learning.rate = 0.001,
                                epsilon.std = 1,
                                null.threshold = 0.05,
                                call.threshold = 0.5,
                                activation.function = "leaky_relu",
                                optimizer = "adam",
                                disable.eager.execution = FALSE,
                                sequence.dictionary = amino.acids[1:20]){
  
  if(length(input.sequences) > number.of.sequences) {
    step <- round(length(input.sequences)/number.of.sequences)
  } else {
    step <- 1
  }
  
  if (disable.eager.execution) {
    tensorflow::tf$compat$v1$disable_eager_execution()
  }
  if(encoder %!in% c("onehotEncoder", "propertyEncoder")) {
    stop("Invalid encoder provided, please select either 'onehotEncoder' or 'propertyEncoder'.")
  }
  
  if(encoder == "propertyEncoder" & !all(sequence.dictionary %in% amino.acids[1:20])){ 
    stop("propertyEncoder method is only available for amino acids, please check the sequence.dictionary.")
  }
  
  if(layers != length(hidden.dims)) {
    stop("Ensure the number of layers matches the vector length of hidden.dims provided.")
  }
  
  optimizer.to.use <- switch(optimizer,
                             "adadelta" = optimizer_adadelta,
                             "adagrad" = optimizer_adagrad,
                             "adam" = optimizer_adam,
                             "adamax" = optimizer_adamax,
                             "ftrl" = optimizer_ftrl,
                             "nadam" = optimizer_nadam,
                             "rmsprop" = optimizer_rmsprop,
                             "sgd" = optimizer_sgd, 
                             stop("Please select a compatible optimizer function in the Keras R implementation."))
  K <- keras::backend()
  
  print("Converting to matrix....")
  if(encoder == "onehotEncoder") {
    sequence.matrix <- onehotEncoder(input.sequences, 
                                     sequence.dictionary = sequence.dictionary,
                                     convert.to.matrix = TRUE)
  } else if (encoder == "propertyEncoder") {
    if(any(aa.method.to.use %!in% names(apex_AA.data))) {
      stop(paste0("Please select one of the following for aa.method.to.use: ", paste(sort(names(apex_AA.data)), collapse = ", ")))
    }
    sequence.matrix <- propertyEncoder(input.sequences, 
                                       method.to.use = aa.method.to.use,
                                       convert.to.matrix = TRUE)
  }
  
  input_shape <- dim(sequence.matrix)[2]
  input_seq <- layer_input(shape = c(input_shape))
  
  # Encoder
  h <- .create_dense_layers(input_seq, 
                            number.of.layers = layers,
                            sizes = hidden.dims,
                            activation.function = activation.function,
                            prefix = "E")
  z_mean <- layer_dense(h, 
                        units = latent.dim, 
                        name = "latent")
  z_log_var <- layer_dense(h, 
                           units = latent.dim, 
                           name = "log_var")
  
  z <- layer_concatenate(list(z_mean, z_log_var)) %>%
       layer_lambda(.vae_sampling,
                    arguments = list(latent.dim = latent.dim, epsilon.std = epsilon.std),
                    name = "lambda")
  
  # Decoder 
  decoder_h <- .create_dense_layers(z, 
                                    number.of.layers = layers,
                                    sizes = rev(hidden.dims),
                                    activation.function = activation.function,
                                    prefix = "D")
  
  decoder_output <- layer_dense(decoder_h, 
                                units = input_shape, 
                                activation = 'sigmoid', 
                                name = "output")
  
  
  # Autoencoder
  vae <- keras_model(inputs = input_seq, 
                     outputs = decoder_output)
  vae %>% keras::compile(optimizer = optimizer.to.use(learning.rate),
                         loss = .vae_loss)
  
  
  print("Fitting Model....")
  # Train the model
  vae %>% fit(
    x = sequence.matrix,
    y = sequence.matrix,
    validation_split = 0.2,
    shuffle = TRUE,
    epochs = epochs,
    batch_size = batch.size,
    #verbose = 0
  )
  
  # Extract Encoder model
  encoder_model <- keras_model(inputs = input_seq, 
                               outputs = list(z_mean, z_log_var))
  
  # Extract Decoder model
  decoder_input <- layer_input(shape = c(latent.dim), name = "decoder_input")
  decoder_h_for_model <- .create_dense_layers(decoder_input, 
                                              number.of.layers = layers,
                                              sizes = rev(hidden.dims),
                                              activation.function = activation.function,
                                              prefix = "D")
  
  decoder_output_for_model <- layer_dense(decoder_h_for_model, 
                                          units = input_shape, 
                                          activation = 'sigmoid', 
                                          name = "output")
  
  decoder_model <- keras_model(inputs = decoder_input, 
                               outputs = decoder_output_for_model)
  
  sequences_encoded <- encoder_model %>% 
                          predict(sequence.matrix, 
                                  batch_size = batch.size)
  
  #TODO allow for variation/n 
  #TODO call for sepcific number of sequences
  z_mean <- K$cast(sequences_encoded[[1]],  "float64")
  z_log_var <- K$cast(sequences_encoded[[2]], "float64")
  
  # Generate epsilon with the same dtype as z_mean and z_log_var
  eps <- k_random_normal(shape = dim(z_mean), 
                         mean = 0, 
                         stddev = 1, 
                         dtype = "float64")
  
  # Sample from the latent space
  z_sample <- z_mean + k_exp(z_log_var / 2) * eps
  
  decoded.sequences <- decoder_model %>% 
    predict(z_sample, 
            steps = step,
            batch_size = batch.size)
  
  new.sequences <- sequenceDecoder(decoded.sequences,
                                   encoder = encoder,
                                   aa.method.to.use = aa.method.to.use,
                                   call.threshold = call.threshold,
                                   sequence.dictionary = sequence.dictionary,
                                   padding.symbol = ".")
  
  return(new.sequences)
}

# Loss and Compilation
#' @importFrom keras k_mean
.reconstruction_loss <- function(y_true, y_pred) {
  k_mean(loss_binary_crossentropy(y_true, y_pred), axis = c(-1))
}

#' @importFrom keras k_square k_exp k_sum
.kl_loss <- function(z_mean, z_log_var) {
  -0.5 * k_sum(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1)
}

#' @importFrom keras k_mean
.vae_loss <- function(y_true, y_pred) {
  reconstruction_loss <- .reconstruction_loss(y_true, y_pred)
  kl_loss <- .kl_loss(y_true, y_pred)
  k_mean(reconstruction_loss + kl_loss)
}

#' @importFrom keras k_random_normal k_shape k_exp
.vae_sampling <- function(arg, latent.dim, epsilon.std){
  z_mean <- arg[, 1:(latent.dim)]
  z_log_var <- arg[, (latent.dim + 1):(2 * latent.dim)]
  
  epsilon <- keras::k_random_normal(
    shape = c(keras::k_shape(z_mean)[[1]]),
    mean=0.,
    stddev=epsilon.std
  )
  
  z_mean + keras::k_exp(z_log_var/2)*epsilon
}

#' @importFrom keras layer_dense
.create_dense_layers <- function(input, 
                                 number.of.layers,
                                 sizes,
                                 activation.function,
                                 prefix) {
  h <- input
  for (i in seq_len(number.of.layers)) {
    h <- layer_dense(h, units = sizes[i], activation = activation.function, name = paste(prefix, i, sep = "."))
  }
  return(h)
}
