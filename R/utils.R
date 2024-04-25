amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", ".")

#Add additional sequence padding to max length
.padded.strings <- function(strings, 
                            max.length,
                            padded.token = NULL,
                            concatenate = TRUE) {
  max_length <- max.length
  
  #Produce combined padded strings
  if(concatenate) {
    x <- lapply(strings, function(str) {
      str_len <- nchar(str)
      if (str_len < max_length) {
        str <- paste(c(str, rep(padded.token, max_length - str_len)), collapse = "")
      } else {
        str
      }
    })
  #Produce a list of strings with seperate vectors
  } else {
    x <- lapply(strings, function(str) {
      str_len <- length(str)
      if (str_len < max_length) {
        str <- c(str, rep(padded.token, max_length - str_len))
      } else {
        str
      }
    })
  }
  return(x)
  
}

#Split Strings by motif.length
substring.extractor <- function(strings, motif.length) {
  lapply(strings, function(x) {
    # Determine the length of the current string
    string_length <- nchar(x)
    
    # Calculate the number of substrings possible
    num_substrings <- string_length - motif.length + 1
    
    if (num_substrings > 0) {
      # Generate all substrings of the specified length
      substrings <- sapply(1:num_substrings, function(j) {
        substring(x, j, j + motif.length - 1)
      })
    } else {
      # Return NA if the string is too short
      substrings <- NA
    }
  }) -> motif.list
  return(motif.list)
}


