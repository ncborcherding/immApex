#' A list of amino acid properties
#' 
#' @description A list of amino acid properties that are
#' used for ```property.encoder()``` function.
#'  
#' @docType data
#' @name apex_AA.data
#' 
NULL

#' A list of IMGT gene names by genes, loci, and species
#' 
#' @description A list of gene names to use for converting for 
#' uniformity. Data pulled from IMGT.org 09/05/2024.
#'  
#' @docType data
#' @name apex_gene.list
#' 
NULL

#' Example contig data for Apex
#' 
#' @description A list of example contigs from the following sources:
#' 
#' \itemize{
#'   \item{TenX: 10k_Human_DTC_Melanoma_5p_nextgem_Multiplex from 10x Website}
#'   \item{AIRR: Human_colon_16S8157851 from \href{https://pubmed.ncbi.nlm.nih.gov/37055623/}{PMID: 37055623}}
#'   \item{Adaptive: Adaptive_2283_D0 from \href{https://pubmed.ncbi.nlm.nih.gov/36220826/}{PMID: 36220826}}
#'   \item{Omniscope: Internal Data}
#'   }
#'  
#' @docType data
#' @name apex_example.data
#' 
NULL

#' List of amino acid substitution matrices
#' 
#' @description A list of amino acid substitution matrices, including:
#' 
#' * BLOSUM45
#' * BLOSUM50  
#' * BLOSUM62  
#' * BLOSUM80  
#' * BLOSUM100
#' * PAM30
#' * PAM40
#' * PAM70    
#' * PAM120
#' * PAM250   
#' 
#' @docType data
#' @name apex_blosum.pam.matrices
#' 