getIMGT <- function(species = "Homo sapiens",
                    chain = "TRB",
                    sequence.type = "aa"
                    region = "v", 
                    release = 7.6) {
  
}

if(region %!in% c("v", "d", "j", "c")) {
  stop("Please select a region in the following category: v, d, j, c")
}



F+ORF+all P = 7.2
F+ORF+in-frame P = 7.5
F+ORF+in-frame P AA = 7.6
F+ORF+in-frame P with IMGT gaps 7.1
F+ORF+in-frame P with IMGT gaps AA 7.3

chain.update <- toupper(paste0("+", chain, region))
release <- 7.6
species.update <- stringr::str_replace(species, " ", "+")

base.url <- https://www.imgt.org/genedb/GENElect?query=7.6+IGHV&species=Homo+sapiens
base.url <- "https://www.imgt.org/genedb/GENElect?query="

updated.url <- paste0(base.url, release, chain.update, "&species=", species.update)


.parse.url <- function(x) {
  if not in_args['gapped'] and not in_args['in_frame_p']:
    return url_stem + "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP", \
  ['IMGT-GENEDB', 'ungapped', 'allP']
  
  elif not in_args['gapped'] and in_args['in_frame_p']:
    return url_stem + "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP", \
  ['IMGT-GENEDB', 'ungapped', 'inframeP']
  
  elif in_args['gapped'] and in_args['in_frame_p']:
    return url_stem + "IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP", \
  ['IMGT-GENEDB', 'gapped', 'inframeP']
  
  elif in_args['gapped'] and not in_args['in_frame_p']:
    raise IOError("IMGT/GENE-DB does not offer a gapped nucleotide FASTA with in frame pseudogenes.")
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

