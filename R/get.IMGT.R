#' Get IMGT Sequences for Specific Loci
#' 
#' Use this to access the ImMunoGeneTics (IMGT) sequences for a
#' specific species and gene loci. More information on 
#' IMGT can be found at {https://www.imgt.org/}{imgt.org}.
#' 
#' @examples
#'TRBV_aa <- get.IMGT(species = "human",
#'                    chain = "TRB",
#'                    frame = "inframe",
#'                    region = "v",
#'                    sequence.type = "aa") 
#'                               
#' @param species One or two word designation of species. Currently supporting: 
#' "human", "mouse", "rat", "rabbit", "rhesus monkey", "sheep", "pig", "platypus",
#' "alpaca", "dog", "chicken", and "ferret"
#' @param chain Sequence chain to access
#' @param frame Designation for "all", "inframe" or "inframe+gap"
#' @param region Sequence gene loci to access
#' @param sequence.type Type of sequence - "aa" for amino acid or "nt" for nucleotide
#' 
#' @importFrom stringr str_extract str_replace_all str_remove_all str_split
#' @importFrom httr GET
#' @importFrom rvest read_html html_text
#' @export get.IMGT
#' @return A list of allele sequences
get.IMGT <- function(species = "human",
                     chain = "TRB",
                     sequence.type = "aa",
                     frame = "inframe",
                     region = "v") {

  if(tolower(region) %!in% c("v", "d", "j", "c")) {
    stop("Please select a region in the following category: 'v', 'd', 'j', 'c'")
  }
  
  if(tolower(sequence.type) %!in% c("aa", "nt")) {
    stop("Please select a sequence.type in the following category: 'aa' or 'nt'")
  }
  
  if(tolower(frame) %!in% c("all", "inframe", "inframe+gap")) {
    stop("Please select a frame in the following category: 'all', 'inframe', or 'inframe+gap'")
  }
  
  if(frame == "inframe+gap" & chain %in% c("IGLJ", "IGKJ", "IGHJ", "IGHD")) {
    stop("IMGT-gapped sequences are not available for 'IGLJ', 'IGKJ', 'IGHJ' or 'IGHD'")
  }
  
  
  selection <- paste0(tolower(frame), "_", tolower(sequence.type))
  
  selection <- switch(selection, 
                     "all_nt" = 7.2,
                     "inframe_nt" = 7.5,
                     "inframe_aa" = 7.6,
                     "inframe+gap_nt" = 7.1,
                     "inframe+gap_aa" = 7.3,
                     "The selection made for sequence.type and frame is not available.")
  
  #Formatting selection into URL for IMGT fasta
  chain.update <- toupper(paste0("+", chain, region))
  species.update <- parseSpecies(species)
  species.update <- stringr::str_replace_all(species.update, " ", "+")
  base.url <- "https://www.imgt.org/genedb/GENElect?query="
  updated.url <- paste0(base.url, selection, chain.update, "&species=", species.update)
  
  print("Getting the sequences from IMGT...")
  response <- httr::GET(updated.url)
  webpage <- content(response, as = "text")
  html <- read_html(webpage)
  pre_text <- html_text(html_nodes(html, "pre"))[2]
  pre_text <- str_remove_all(pre_text, "\n")
  sequences <- str_split(pre_text, ">")[[1]]
  sequences <- sequences[nchar(sequences) > 0]  # remove any empty entries
  
  print("Formatting IMGT sequences...")
  fasta_list <- list()
  for (seq in sequences) {
    # Extract the allele name using regex. Allele name generally follows the pattern before the first "|"
    parts <- str_split(seq, "\\|")[[1]]
    if (length(parts) >= 2) {
      allele_name <- parts[2]  # This should correspond to something like 'TRBV1*01'
    } else {
      next  # Skip this sequence if it does not have enough parts
    }
    
    # Remove any non-sequence characters (like digits, description text, etc.)
    if(tolower(sequence.type) == "aa") {
      sequence <- str_extract(seq, "(?<=\\|)[\\s\\S]*$")
      sequence <- gsub("[^A-Z]", "", sequence)  # Assuming only uppercase letters are valid
      sequence <- stringr::str_remove_all(sequence, paste0(chain, toupper(region)))
    } else if(tolower(sequence.type) == "nt") {
      sequence <- str_extract(seq, "[acgt]+$")
      sequence <- gsub("[^a-z]", "", sequence)  # Assuming only uppercase letters are valid
    }
    
    # Assign the sequence to the allele name in the list
    fasta_list[[allele_name]] <- sequence
  }
  
  result.list <- list(sequences = fasta_list,
                      misc = list(species = species.update,
                                  chain = chain,
                                  sequence.type = sequence.type,
                                  frame = frame,
                                  region = region))
  return(result.list)

}

#' @importFrom hash hash
parseSpecies <- function(x) {
  
  species <- c("human", "mouse", "rat", "rabbit",
               "rhesus monkey", "sheep", "pig", "platypus",
               "alpaca", "dog", "chicken", "ferret")
  
  species_dictionary <- hash::hash(
    "human" = "Homo sapiens",
    "mouse" = "Mus",
    "rat" = "Rattus norvegicus",
    "rabbit" = "Oryctolagus cuniculus",
    "rhesus monkey" = "Macaca mulatta",
    "sheep" = "Ovis aries",
    "pig" = "Sus scrofa",
    "platypus" = "Ornithorhynchus anatinus",
    "alpaca" = "Vicugna pacos",
    "dog" = "Canis lupus familiaris",
    "chicken" = "Gallus gallus",
    "ferret" = "Mustela putorius furo"
  )
  
  x <- tolower(x)
  if(x %in% species) {
    return(species_dictionary[[x]])
  } else {
    stop(paste0("Please select one of the following species: ", paste(species, collapse = ", ")))
  }
}

