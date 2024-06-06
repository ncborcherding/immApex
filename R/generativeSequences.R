generativeSequences <-function(input.sequences,
                               encoder = "onehotEncoder",
                               aa.method.to.use = NULL,
                               layers = 2, 
                               hidden.dims = c(256, 128),
                               latent.dim = 16,
                               batch.size = 16,
                               epochs = 30,
                               learning.rate = 0.0001,
                               epsilon.std = 1,
                               null.threshold = 0.05,
                               activation.function = "relu",
                               optimizer = "adam",
                               sequence.dictionary = amino.acids[1:20]){
 
  
  if (tensorflow::tf$executing_eagerly()) {
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
  
  if(encoder == "onehotEncoder") {
    sequence.matrix <- onehotEncoder(input.sequences, 
                                     sequence.dictionary = sequence.dictionary,
                                     convert.to.matrix = TRUE)
  } else if (encoder == "propertyEncoder") {
    if(any(method.to.use %!in% names(apex_AA.data))) {
      stop(paste0("Please select one of the following for aa.method.to.use: ", paste(sort(names(apex_AA.data)), collapse = ", ")))
    }
    sequence.matrix <- propertyEncoder(input.sequences, 
                                       sequence.dictionary = aa.method.to.use,
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
  z_mean <- layer_dense(h, units = latent.dim, 
                        name = "latent")
  z_log_var <- layer_dense(h, units = latent.dim, 
                           name = "log_var")
  
  z <- layer_concatenate(list(z_mean, z_log_var)) %>% 
       layer_lambda(.vae_sampling, name = "lambda")
  
  #Decoder 
  decoder_h <- .create_dense_layers(z, 
                                    number.of.layers = layers,
                                    sizes = rev(hidden.dims),
                                    activation.function = activation.function,
                                    prefix = "D")
  
  decoder_output <- layer_dense(decoder_h, 
                                units = input_shape, 
                                activation = 'sigmoid', 
                                name = "output")
  #Autoencoder
  vae <- keras_model(inputs = input_seq, 
                     outputs = decoder_output)
  vae %>% keras::compile(optimizer = optimizer.to.use(learning.rate),
                         loss = .vae_loss)
  
  
  # Train the model
  vae %>% fit(
    x = sequence.matrix,
    y = sequence.matrix,
    validation_split = 0.2,
    shuffle = TRUE,
    epochs = epochs,
    batch_size = batch.size
  )
}
  
output <- layer_dense(dec, units = 1, activation = 'sigmoid', name = "output")

# Create the model
model <- keras_model(inputs = decoder_input, outputs = decoder_output)

# Loss and Compilation
.reconstruction_loss <- function(y_true, y_pred) {
  k_mean(loss_binary_crossentropy(y_true, y_pred), axis = c(-1))
}

.kl_loss <- function(y_true, y_pred) {
  -0.5 * k_sum(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1)
}

.vae_loss <- function(y_true, y_pred) {
  .reconstruction_loss(y_true, y_pred) + .kl_loss(y_true, y_pred)
}

.vae_sampling <- function(arg){
  z_mean <- arg[, 1:(latent.dim)]
  z_log_var <- arg[, (latent.dim + 1):(2 * latent.dim)]
  
  epsilon <- keras::k_random_normal(
    shape = c(keras::k_shape(z_mean)[[1]]),
    mean=0.,
    stddev=epsilon.std
  )
  
  z_mean + keras::k_exp(z_log_var/2)*epsilon
}

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
