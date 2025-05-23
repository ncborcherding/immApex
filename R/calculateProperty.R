#' Position-wise Amino-Acid Property Profiles
#'
#' Computes a range of summary statistics for property values of one or more AA 
#' property scales at ever esidue position of a set of protein (or peptide) 
#' sequences.  The function is entirely vectorised: it first calls 
#' [`calculateFrequency()`] to obtain a residue-by-position **frequency** 
#' matrix *F* (each column sums to 1) and then performs a single matrix product.
#'
#' @param sequences Character vector of amino-acid strings.
#' @param property.set See [`positionalPropertyProfile()`].
#' @param summary.fun Character string (`"mean"`, `"median"`, `"sum"`,
#' `"min"`, `"max"`), **or** a function accepting a numeric vector and
#' returning length-1 numeric.  Defaults to `"mean"`.
#' @param transform Character string controlling a *post-summary*
#' transformation.  One of  `"none"` (default), `"sqrt"`, `"log1p"`, 
#' `"zscore"` (row-wise), or `"minmax"` (row-wise).
#' @param groups Optional group factor (see [`calculateProperty()`]).  
#' Summary statistics are computed independently within each group.
#' @param max.length Integer.  Pad/trim to this length
#'   (`max(nchar(sequences))` by default).
#' @param padding.symbol Single character used for right-padding. Must not be
#'   one of the 20 canonical residues.
#' @param tidy Logical; if `TRUE`, return a long-format `data.frame`
#'
#' @return
#' *When `groups` is **NULL*** – a numeric matrix (*k* × *L*).  
#' *When `groups` is supplied* – an array of dimension
#'   (*k* × *L* × *G*) with dimnames `scale`, `position`, `group`.  
#' If `tidy = TRUE`, a long `data.frame` is returned instead
#' .
#' @export
positionalPropertySummary <- function(sequences,
                                      property.set = "Atchley",
                                      summary.fun  = "mean",
                                      transform    = "none",
                                      groups       = NULL,
                                      max.length   = NULL,
                                      padding.symbol = ".",
                                      tidy           = FALSE) {
  
  ## ------------------------------------------------------------------ 
  ## 0.  Helpers + argument normalization 
  ## ------------------------------------------------------------------ 
  stopifnot(is.character(sequences),
            padding.symbol %!in% amino.acids,
            nchar(padding.symbol) == 1L)
  
  # coerce summary.fun into function (fast path for built-ins)
  summary.fun <- switch(
    deparse(substitute(summary.fun)),
    "mean"  = base::mean,
    "median"= stats::median,
    "sum"   = base::sum,
    "min"   = base::min,
    "max"   = base::max,
    summary.fun)   # user-supplied
  
  if (!is.function(summary.fun) || length(formals(summary.fun)) != 1)
    stop("'summary.fun' must be a one-argument function or one of the built-in names.")
  
  transform <- match.arg(transform,
                         c("none","sqrt","log1p","zscore","minmax"))
  
  if (!is.null(groups)) {
    groups <- as.factor(groups)
    stopifnot(length(groups) == length(sequences))
  }
  
  if (is.null(max.length))
    max.length <- max(nchar(sequences), 1L)
  
  ## ------------------------------------------------------------------ ##
  ## 1.  Property matrix (k × 20) ------------------------------------ ##
  ## ------------------------------------------------------------------ ##
  # -- (a) built-in scale sets
  if (is.character(property.set)) {
    
    S <- as.matrix(S)
    
    # -- (b) user-supplied matrix ----------------------------------------
  } else if (is.matrix(property.set) || is.data.frame(property.set)) {
    
    S <- as.matrix(property.set)
    
  } else {
    stop("'property.set' must be a recognised name or a numeric matrix.")
  }
  
  # -- (c) sanity checks & AA column ordering ---------------------------
  if (is.null(colnames(S)))
    stop("Custom property matrices must have the 20 amino acids as column names.")
  
  missing <- setdiff(amino.acids, colnames(S))
  if (length(missing))
    stop("Property matrix is missing columns for: ",
         paste(missing, collapse = ", "))
  
  # TODO: write .aa_property_matrix
  #S <- .aa_property_matrix(property.set)
  S <- switch(property.set,
              "Atchley"  = protr::extractAtchleyFactor(),
              "Kidera"   = protr::extractKideraFactor(),
              "stScales" = protr::extractStScales(),
              "tScales"  = protr::extractTScales(),
              "VHSE"     = protr::extractVHSEScales())
  
  S <- S[ , amino.acids, drop = FALSE]   
  k <- nrow(S)                            
  
  ## ------------------------------------------------------------------ 
  ## 2.  Generate residue counts
  ## ------------------------------------------------------------------
  counts_per_group <- function(seq_subset) {
    pad_mat <- .padded.strings(seq_subset,
                               max.length = max.length,
                               padded.token = padding.symbol,
                               concatenate = TRUE)
    seq_mat <- do.call(rbind, strsplit(unlist(pad_mat), ""))
    
    # Map to 1:21 integers 
    lvl <- c(amino.acids, padding.symbol)
    idx <- match(seq_mat, lvl)
    
    # Per-position tabulation (20 rows × L cols) 
    C <- matrix(0L, 20, max.length,
                dimnames = list(amino.acids, NULL))
    for (j in seq_len(max.length))
      C[ , j] <- tabulate(idx[ , j], nbins = 21L)[1:20]
    C
  }
  
  grp_levels <- if (is.null(groups)) "(all)" else levels(groups)
  Summ <- array(NA_real_, dim = c(k, max.length, length(grp_levels)),
                dimnames = list(scale    = rownames(S),
                                position = paste0("Pos.", seq_len(max.length)),
                                group    = grp_levels))
  
  for (g in seq_along(grp_levels)) {
    
    idx_g <- if (is.null(groups)) seq_along(sequences) else which(groups == grp_levels[g])
    C     <- counts_per_group(sequences[idx_g])         
    Npos  <- colSums(C)
    
    if (identical(summary.fun, base::mean)) {          
      Freq <- sweep(C, 2L, Npos, "/")
      Summ[ , , g] <- S %*% Freq                         
      next
    }
    
    if (identical(summary.fun, base::sum)) {             
      Summ[ , , g] <- S %*% C                         
      next
    }
    
    ## ---------------------------------------------------------------------- 
    ## 3.  Generic summaries (median, min, max, custom) 
    ## ---------------------------------------------------------------------- 
    # Pre-compute, for each scale, the property value per residue
    prop_by_res <- lapply(seq_len(k), function(i) S[i, ])
    
    for (j in seq_len(max.length)) {
      cnts <- C[ , j]                                    # length-20
      if (!any(cnts)) next                               # all padding
      for (i in seq_len(k)) {
        vals <- rep(prop_by_res[[i]], times = cnts)
        Summ[i, j, g] <- summary.fun(vals)
      }
    }
  }
  
  ## ------------------------------------------------------------------ 
  ## 4.  Optional row-wise transformation 
  ## ------------------------------------------------------------------ 
  transform_apply <- switch(transform,
                            "none"   = identity,
                            "sqrt"   = function(x) sqrt(pmax(x, 0)),
                            "log1p"  = function(x) log1p(pmax(x, 0)),
                            "zscore" = function(x) {
                              center <- rowMeans(x, na.rm = TRUE)
                              scale  <- sqrt(pmax(rowMeans(x^2, na.rm = TRUE) - center^2, .Machine$double.eps))
                              sweep(sweep(x, 1L, center), 1L, scale, "/")
                            },
                            "minmax" = function(x) {
                              a <- apply(x, 1L, min,  na.rm = TRUE)
                              b <- apply(x, 1L, max,  na.rm = TRUE)
                              sweep(sweep(x, 1L, a), 1L, pmax(b - a, .Machine$double.eps), "/")
                            })
  Summ <- transform_apply(Summ)
  
  ## ------------------------------------------------------------------ 
  ## 5.  Tidy reshaping  
  ## ------------------------------------------------------------------ 
  if (tidy) {
    if (length(grp_levels) == 1L) {
      Summ <- as.data.frame(as.table(Summ[ , , 1, drop = FALSE]),
                            stringsAsFactors = FALSE,
                            responseName = "value")
      names(Summ) <- c("scale","position","value")
    } else {
      Summ <- reshape2::melt(Summ, varnames = c("scale","position","group"),
                             value.name = "value")
    }
    Summ$position <- as.integer(sub("Pos\\.", "", Summ$position))
  }
  
  Summ
}
