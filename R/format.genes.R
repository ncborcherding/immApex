#' @importFrom stringr str_split
format.genes <- function(input.data,
                         genes = "v",
                         technology = NULL,
                         species = "human",
                         simplify.format = TRUE) {
  
  if(technology %in% c("TenX","Adaptive")) {
    genes.updated <- paste0(genes, "_gene")
    if(any(genes.updated %!in% colnames(input.data))) {
      genes.updated <- paste0(genes, "MaxResolved")
    }
  } else if (technology %in% c("AIRR", "Omniscope")) {
    genes.updated <- paste0(genes, "_call")
  }
  #TODO Add screpertoire support
  
  gene.list <- apex_gene.list[tolower(genes)]
  lapply(gene.list, function(x) {
    lapply(x, function(y) {
      reference <- y[[tolower(species)]]
      if(simplify.format) {
        reference <- unique(str_split(reference, "[*]", simplify = TRUE)[,1])
      }
      reference
    }) -> segment.reference
  }) -> gene.reference
  
  for(i in seq_along(genes)) {
    input.data[,genes.updated[i]] <- str_split(input.data[,genes.updated[i]], "[|]", simplify = TRUE)[,1]
    if(technology == "Adaptive") {
      input.data[,genes.updated[i]] <- .remove_leading_zeros(input.data[,genes.updated[i]], genes[i])
    }
    if(simplify.format) {
      input.data[,genes.updated[i]] <- str_split(input.data[,genes.updated[i]], "[*]", simplify = TRUE)[,1]
    }
    #TODO Need to method match simplified format or more complex
    #TODO Test AIRR, Omniscope, 10x
    which(input.data[,genes.updated[i]] %!in% unlist(gene.reference[[i]]))
    
  
  }

#Working on adaptive gene formatig
#' @importFrom stringr str_replace_all
.remove_leading_zeros <- function(x, 
                                  gene) {
    #Mod Gene name to remove TCR
    x <- str_replace_all(x, "TCR", "TR") 
    pattern <- toupper(paste0(gene, "0"))
    
    # Remove leading zeros after gene
    x <- str_replace_all(x, pattern, toupper(gene))
    
    # Split at "-", remove leading zeros from the second part, and paste back
    sapply(strsplit(x, "-"), function(parts) {
      if (length(parts) == 2) {
        paste0(parts[1], "-", sub("^0+", "", parts[2]))
      } else {
        parts[1]
      }
    })
}
  
 

reference <- as.data.frame(readxl::read_xlsx("~/Documents/GitHub/OS-toolbox/data/gene.reference2.xlsx"))
