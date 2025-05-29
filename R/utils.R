amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

"%!in%" <- Negate("%in%")

.check.sequences <- function(sequences, dict) {
  pat <- sprintf("[^%s]", paste(dict, collapse = ""))
  any(grepl(pat, head(sequences, 20L)))  
}

#Add additional sequence padding to max length
.padded.strings <- function(strings, max.length, pad = "", collapse = TRUE) {
  
  pad_right  <- sprintf("%%-%ds", max.length)          # e.g. "%-20s"
  pad_vector <- function(s) c(strsplit(s, "")[[1L]],
                              rep(pad, max.length - nchar(s)))
  
  if (collapse) {
    sprintf(pad_right, strings)                        # base-R, no loop
  } else {
    lapply(strings, function(s) pad_vector(s))
  }
}

.substring.extractor <- function(strings, k) {
  vapply(strings, function(s) {
    n <- nchar(s)
    if (n < k) return(NA_character_)
    starts <- seq_len(n - k + 1L)
    substring(s, starts, starts + k - 1L)
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
}


.min.max.normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))   # avoids div-by-zero
  (x - rng[1L]) / diff(rng)
}


#' @importFrom matrixStats colMedians colMeans2 colSums2 colVars colMads
.get.stat.function <- function(method) {
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("matrixStats is required for summary statistics.", call. = FALSE)
  switch(method,
         median = matrixStats::colMedians,
         mean   = matrixStats::colMeans2,
         sum    = matrixStats::colSums2,
         vars   = matrixStats::colVars,
         mads   = matrixStats::colMads,
         stop("Invalid `method`"))
}

.is_seurat_object <- function(obj) inherits(obj, "Seurat")
.is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
.is_seurat_or_se_object <- function(obj) {
  .is_seurat_object(obj) || .is_se_object(obj)
}

.get.genes.updated <- function(data, tech, region) {
  if (is_container(data)) return(region)
  
  tech <- tech %||% "other"
  key  <- list(TenX = paste0(region, "_gene"),
               Adaptive = paste0(region, "_gene"),
               AIRR = paste0(region, "_call"))[[tech]]
  
  if (is.null(key) || !(key %in% colnames(data)))
    key <- paste0(region, "GeneName")
  
  key
}

.array.dimnamer <- function(arr) {
  comb <- expand.grid(dimnames(arr)[2:3], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  comb <- comb[order(as.integer(comb[,1])), ]    # numeric-aware sort
  paste(comb[,1], comb[,2], sep = "_")
}

.array.reshape <- function(x, dim, order = c("C", "F")) {
  
  order <- match.arg(order)
  dim   <- as.integer(dim)
  
  # Checks 
  if (!is.atomic(x))
    stop("`x` must be an atomic vector, matrix, or array.", call. = FALSE)
  if (any(dim <= 0))
    stop("All elements of `dim` must be positive.", call. = FALSE)
  if (length(x) != prod(dim))
    stop(sprintf("Input of length %d cannot be reshaped to %s (product = %d).",
                 length(x), paste(dim, collapse = "x"), prod(dim)),
         call. = FALSE)
  
  # row-major flatten 
  row_flatten <- function(a) {
    if (is.null(dim(a))) {
      as.vector(a)                                 # already 1D
    } else {
      as.vector(aperm(a, rev(seq_along(dim(a)))))  # reverse dims then flatten
    }
  }
  
  #  build result 
  if (order == "F") {
    array(as.vector(x), dim = dim)
    
  } else {  
    vals_C <- row_flatten(x)
    tmp    <- array(vals_C, dim = rev(dim))
    aperm(tmp, rev(seq_along(dim)))
  }
}

.get_v_column <- function(input.data, technology = c("TenX", "Adaptive", "AIRR")) {
  stopifnot(is.data.frame(input.data))
  technology <- match.arg(technology)
  
  if ("v_IMGT" %in% names(input.data)) {
    v.col <- "v_IMGT"
  } else {
    v.col <- switch(technology,
                    TenX     = if ("v_gene" %in% names(input.data)) "v_gene" else "vGeneName",
                    Adaptive = if ("v_gene" %in% names(input.data)) "v_gene" else "vGeneName",
                    AIRR     = "v_call")
    )
  }
  
  if (!v.col %in% names(input.data)) {
    stop("Cannot find a V-gene column in `input.data` for technology = '", technology, "'.")
  }
  
  v.col
}