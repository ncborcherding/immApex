#' Ensure clean gene nomenclature using IMGT annotations
#' 
#' This function will format the genes into a clean
#' nomenclature using the IMGT conventions. 

#' @param input.data Data frame of sequencing data or 
#' scRepertoire outputs
#' @param region Sequence gene loci to access - 'v', 'd', 'j', 'c'
#' or a combination using c('v', 'd', 'j')
#' @param technology The sequencing technology employed - 'TenX', "Adaptive', 
#' 'AIRR', or 'Omniscope'.
#' @param species One or two word designation of species. Currently supporting: 
#' "human", "mouse", "rat", "rabbit", "rhesus monkey", "sheep", "pig", "platypus",
#' "alpaca", "dog", "chicken", and "ferret"
#' @param simplify.format If applicable, remove the allelic designation (TRUE) or
#' retain all information (FALSE)
#' 
#' @examples
#' format.genes(apex_example.data[["TenX"]],
#'              region = "v",
#'              technology = "TenX")
#' 
#' @importFrom stringr str_split
#' 
#' @export format.genes
#' @return A data frame with the new columns of formatted genes added.


#' @importFrom stringr str_split
format.genes <- function(input.data,
                         region = "v",
                         technology = NULL,
                         species = "human",
                         simplify.format = TRUE) {
  data("apex_gene.list")
  if(any(tolower(region) %!in% c("v", "d", "j", "c"))) {
    stop("Please select a region or regions in the following category: 'v', 'd', 'j', 'c'")
  }
  if(!.is_seurat_or_se_object(input.data)) {
    if(technology %!in% c("TenX", "AIRR", "Adaptive", "Omniscope")) {
      stop("Please select a technology in the following category: 'TenX', 'AIRR', 'Adaptive', 'Omniscope'")
    }
  }
  
  if (.is_seurat_or_se_object(input.data)) {
    chain.1 <- getIR(input.data, 
                     chains = "TRA", 
                     sequence.type = "aa")[[1]]
    chain.2 <- getIR(input.data, 
                     chains = "TRB", 
                     sequence.type = "aa")[[1]]
    input.data <- rbind(chain.1, chain.2)
    genes.updated <- region
  } else {
    input.data[input.data == ""] <- NA
    if(technology %in% c("TenX","Adaptive")) {
      genes.updated <- paste0(region, "_gene")
      if(any(genes.updated %!in% colnames(input.data))) {
        genes.updated <- paste0(region, "GeneName")
        input.data[,genes.updated][is.na(input.data[,genes.updated])] <- str_split(input.data[,"vGeneNameTies"][is.na(input.data[,genes.updated])], "[,]", simplify = TRUE)[,1]
      }
    } else if (technology %in% c("AIRR", "Omniscope")) {
      genes.updated <- paste0(region, "_call")
    }
  }
 
  
  gene.list <- apex_gene.list[tolower(region)]
  lapply(gene.list, function(x) {
    lapply(x, function(y) {
      reference <- y[[tolower(species)]]
      if(simplify.format) {
        reference <- unique(str_split(reference, "[*]", simplify = TRUE)[,1])
      }
      reference
    }) -> segment.reference
  }) -> gene.reference
  
  for(i in seq_along(region)) {
    input.data[,paste0(region[i], "_IMGT")] <- str_split(input.data[,genes.updated[i]], "[|]", simplify = TRUE)[,1]
    if(technology == "Adaptive") {
      input.data[,paste0(region[i], "_IMGT")] <- .remove_leading_zeros(input.data[,paste0(region[i], "_IMGT")], region[i])
    }
    if(simplify.format) {
      input.data[,paste0(region[i], "_IMGT")]  <- str_split(input.data[,paste0(region[i], "_IMGT")], "[*]", simplify = TRUE)[,1]
    }
    input.data[,paste0(region[i], "_IMGT.check")] <- 1
    input.data[,paste0(region[i], "_IMGT.check")][which(input.data[,paste0(region[i], "_IMGT")] %!in% unlist(gene.reference[[i]]))] <- 0
  }
  return(input.data)
}

#Working on adaptive gene formatting
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
