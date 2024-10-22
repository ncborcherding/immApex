#' Infer portions of the CDR Loop from Vgene data
#' 
#' Use this isolate sequences from the CDR loop using
#' the V gene annotation. When there are multiple 
#' V gene matches for a single gene, the first allelic 
#' sequence is used.
#' 
#' @examples
#' # Getting the Sequence Reference
#' data(immapex_example.data)
#' TRBV_aa <- getIMGT(species = "human",
#'                    chain = "TRB",
#'                    frame = "inframe",
#'                    region = "v",
#'                    sequence.type = "aa") 
#'                     
#' # Ensuring sequences are formatted to IMGT                   
#' TenX_formatted <- formatGenes(immapex_example.data[["TenX"]],
#'                               region = "v",
#'                               technology = "TenX")
#'              
#' # Inferring CDR loop elements           
#' TenX_formatted <- inferCDR(TenX_formatted,
#'                            chain = "TRB", 
#'                            reference = TRBV_aa,
#'                            technology = "TenX", 
#'                            sequence.type = "aa",
#'                            sequences = c("CDR1", "CDR2"))
#'                             
#' @param input.data Data frame output of \code{\link{formatGenes}}
#' @param reference IMGT reference sequences from \code{\link{getIMGT}}
#' @param chain Sequence chain to access, like \strong{TRB} or \strong{IGH}
#' @param technology The sequencing technology employed - \strong{TenX}, \strong{Adaptive}, 
#' \strong{AIRR}, or \strong{Omniscope}
#' @param sequence.type Type of sequence - \strong{aa} for amino acid or \strong{nt} for nucleotide
#' @param sequences The specific regions of the CDR loop to get from the data, such as \strong{CDR1}. 
#'                             
#' @importFrom stringr str_split
#' 
#' @export
#' @return A data frame with the new columns of CDR sequences added.
inferCDR <- function(input.data, 
                     reference = NULL,
                     chain = "TRB", 
                     technology = NULL, 
                     sequence.type = "aa",
                     sequences = c("CDR1", "CDR2")) {
  
  region <- "v"
  if(is.null(reference) || "v" %!in% reference[["misc"]][["region"]]) {
    stop("Please provide a list of V gene reference sequence using 'get.IMGT()'.")
  }
  
  if(reference[["misc"]][["sequence.type"]] != sequence.type) {
    stop("Please check if the reference provided matches the sequence.type selected.")
  }
  
  if(technology %!in% c("TenX", "AIRR", "Adaptive", "Omniscope")) {
      stop("Please select a technology in the following category: 'TenX', 'AIRR', 'Adaptive', 'Omniscope'")
  }
  
  if ("v_IMGT" %!in% colnames(input.data)) {
    warning("No output from 'format.genes()' detecting, proceeding with gene nomenclature that may result in NAs.")
    if(technology %in% c("TenX","Adaptive")) {
      genes.updated <- paste0(region, "_gene")
      if(any(genes.updated %!in% colnames(input.data))) {
        genes.updated <- paste0(region, "GeneName")
        input.data[,genes.updated][is.na(input.data[,genes.updated])] <- str_split(input.data[,"vGeneNameTies"][is.na(input.data[,genes.updated])], "[,]", simplify = TRUE)[,1]
      }
    } else if (technology %in% c("AIRR", "Omniscope")) {
      genes.updated <- paste0(region, "_call")
    }
  } else {
    v.col <- "v_IMGT"
  }
  
  sequence.pos <- .sequence.positions[grep(paste0(sequences, collapse = "|"), names(.sequence.positions))]
  if(sequence.type == "nt") {
    lapply(sequence.pos, function(x) {
        start <- min(x*3)
        end <- max(x*3) + 3
        seq(start, end)
    }) -> sequence.pos
  }
  
  v.genes <- unique(input.data[,v.col])
  lapply(seq_len(length(v.genes)), function(x) {
    ref.pos <- which(names(reference[["sequences"]]) == v.genes[x])
    if(length(ref.pos) == 0) {
      ref.pos <- grep(v.genes[x], names(reference[["sequences"]]))[1]
    }
    if(is.na(ref.pos)) {
      return(c(v.genes[x], rep(NA, length(sequences))))
    } else {
      lapply(sequence.pos, function(y) {
         substring(reference[["sequences"]][[ref.pos]], min(y), max(y))
      }) -> string.list
      return(c(v.genes[x], unlist(string.list)))
    }
  }) -> cdr.sequences.list
  
  cdr.matrix <- do.call(rbind, cdr.sequences.list)
  colnames(cdr.matrix) <- c("v_IMGT", paste0(sequences, "_IMGT"))
  
  input.data <- merge(input.data, cdr.matrix, by.x = v.col, by.y = 1)
  return(input.data)
}

#CDR3 AA sequences by IMGT
.sequence.positions <- list(
  FR1.A = c(1:15), 
  FR1.B = c(16:26),
  CDR1 = c(27:38), 
  FR2.C = c(39:46), 
  FR2.Cp = c(47:55),
  CDR2 = c(56:65), 
  FR3.C = c(66:74), 
  FR3.D = c(75:84),
  FR3.E = c(85:96), 
  FR3.F = c(97:104)
)
