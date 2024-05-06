getIMGT <- function(species = "human",
                    chain = "TRB",
                    sequence.type = "aa"
                    frame = "inframe",
                    region = "v") {

if(region %!in% c("v", "d", "j", "c")) {
  stop("Please select a region in the following category: v, d, j, c")
}

selection <- paste0(frame, "_", sequence.type)

selection <- switch(selection, 
                   "all_nt" = 7.2,
                   "inframe_nt" = 7.5,
                   "inframe_aa" = 7.6,
                   "gap_nt" = 7.1,
                   "gap_aa" = 7.3,
                   "The selection made for sequence.type and frame is not available.")


chain.update <- toupper(paste0("+", chain, region))
species.update <- parseSpecies(species)
species.update <- stringr::str_replace(species.update, " ", "+")
base.url <- "https://www.imgt.org/genedb/GENElect?query="

updated.url <- paste0(base.url, selection, chain.update, "&species=", species.update)

#TODO FASTA download
#TODO FASTA conversion
#TODO Add detection

}

#' @importFrom hash hash
parseSpecies <- function(x) {
  
  species <- c("human", "mouse", "rat", "rabbit",
               "rhesus monkey", "sheep", "pig", "platypus",
               "alpaca", "dog", "chicken", "ferret")
  
  species_dictionary <- hash::hash(
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus",
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

