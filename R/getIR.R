#' Extract Immune Receptor Sequences
#' 
#' Use this to extract immune receptor sequences from a Single-Cell 
#' Object or the output of \link[scRepertoire]{combineTCR} and 
#' \link[scRepertoire]{combineBCR}.

#' @param input.data Single-cell object or the output of 
#' \link[scRepertoire]{combineTCR} and \link[scRepertoire]{combineBCR} from
#' scRepertoire.
#' @param chains Immune Receptor chain to use - \strong{TRA}, 
#' \strong{TRB}, \strong{Heavy}, \strong{Light}
#' @param sequence.type Extract amino acid (\strong{aa}) or 
#' nucleotide (\strong{nt})
#' 
#' @importFrom stringr str_split
#' 
#' @export
#' @return A data frame of nucleotide or amino acid sequences 
getIR <- function(input.data, 
                  chains,
                  sequence.type = "aa") {
  
  col.pos <- switch(sequence.type,
                    "aa" = "CTaa",
                    "nt" = "CTnt"
                    stop("Please select either 'aa' or 'nt' for sequence.type."))
  
  if (inherits(x=input.data, what ="Seurat") | inherits(x=input.data, what ="SingleCellExperiment")) {
    meta <- .grabMeta(input.data)
  } else {
    meta <- do.call(rbind,input.data)
    rownames(meta) <- meta[,"barcode"]
  }
  if(chains %!in% c("TRA", "TRB", "Heavy", "Light")) {
    stop("Please select one of the following chains: 'TRA', 'TRB', 'Heavy', 'Light'")
  }
  tmp <- data.frame(barcode = rownames(meta), 
                    str_split(meta[,"CTaa"], "_", simplify = TRUE), 
                    str_split(meta[,"CTgene"], "_", simplify = TRUE))
  if (length(chains) == 1 && chains != "both") {
    if (chains %in% c("TRA", "IGH")) { #here
      pos <- list(c(2,4))
    } else if (chains %in% c("TRB", "IGL")) { #here
      pos <- list(c(3,5))
    }
  } 
  
  IR <- NULL
  for (i in seq_along(pos)) {
    sub <- as.data.frame(tmp[,c(1,pos[[i]])])
    colnames(sub) <- c("barcode", "cdr3_aa", "genes")
    sub$v <- str_split(sub$genes, "[.]", simplify = TRUE)[,1]
    sub$j <- str_split(sub$genes, "[.]", simplify = TRUE)[,2]
    sub[sub == ""] <- NA
    IR[[i]] <- sub
    sub <- NULL
  }
  names(IR) <- chains
  return(IR)
}

#' @importFrom SingleCellExperiment colData 
.grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  } else if (inherits(x=sc, what ="SingleCellExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[clu] <- "cluster.active.idents"
    } else {
      colnames(meta)[clu] <- "cluster"
    }
  }
  return(meta)
}