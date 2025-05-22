#' Shannon Diversity Index (Entropy)
#'
#' Calculates Shannon’s information entropy (often denoted *H′*) for a set
#' of clone or sequence counts.
#'
#' \deqn{H' \;=\; -\sum_{i = 1}^{S} p_i \ln p_i}
#'
#' where *p*<sub>*i*</sub> = *n*<sub>*i*</sub> / *N* are the relative
#' frequencies (proportions) of each of the *S* distinct categories.
#'
#' @param cnt Numeric vector of non-negative counts (one entry per clone/
#'   residue/OTU).  Zero counts are ignored.
#'
#' @return A single numeric value (≥ 0).  When `cnt` contains exactly one
#'   positive entry the function returns `0`.
#'
#' @seealso [norm_entropy()], [inv_simpson()]
#' @examples
#' counts <- c(A = 12, B = 4, C = 4)
#' shannon_entropy(counts)
#' @export
shannon_entropy <- function(cnt) {
  cnt <- cnt[cnt > 0]
  p   <- cnt / sum(cnt)
  -sum(p * log(p))
}

#' Inverse Simpson Diversity
#'
#' Computes the inverse of Simpson’s concentration index, sometimes written
#' as *1/D*.  This metric emphasises dominant categories.
#'
#' \deqn{1/D \;=\; \frac{1}{\sum_{i} p_i^{\,2}}}
#'
#' @inheritParams shannon_entropy
#' @return Numeric value ≥ 1.  Equals 1 when all observations belong to a
#'   single category.
#' @examples
#' inv_simpson(c(10, 5, 1))
#' @export
inv_simpson <- function(cnt) {
  cnt <- cnt[cnt > 0]
  p   <- cnt / sum(cnt)
  1 / sum(p * p)
}

#' Gini–Simpson Diversity
#'
#' Computes the complement of Simpson’s index (also called the
#' Gini–Simpson index or probability of interspecific encounter):
#'
#' \deqn{1 - \lambda = 1 - \sum_{i} p_i^{\,2}}
#'
#' @inheritParams shannon_entropy
#' @return Value in the interval [0, 1].  Higher numbers indicate greater
#'   heterogeneity.
#' @examples
#' gini_simpson(c(10, 5, 5))
#' @export
gini_simpson <- function(cnt) 1 - 1 / inv_simpson(cnt)

#' Normalised Shannon Entropy
#'
#' Shannon entropy scaled to the interval [0, 1] by its maximum possible
#' value given *S* observed categories:
#'
#' \deqn{H^* = \frac{H'}{\ln S}}
#'
#' (also known as “Shannon evenness”).
#'
#' @inheritParams shannon_entropy
#' @return Numeric value in [0, 1]; `0` when all observations are in a
#'   single category.
#' @examples
#' norm_entropy(c(40, 10, 10, 10))
#' @export
norm_entropy <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S   <- length(cnt)
  if (S == 1L) 0 else shannon_entropy(cnt) / log(S)
}

#' Pielou’s Evenness
#'
#' Convenience wrapper for normalized Shannon entropy (*E* = *H′* / ln *S*).
#'
#' @inheritParams shannon_entropy
#' @return Numeric evenness measure in [0, 1].
#' @examples
#' pielou_evenness(c(3, 3, 3))
#' @export
pielou_evenness <- function(cnt) norm_entropy(cnt)

#' Hill-Number Generator
#'
#' Returns a *function* that computes the Hill diversity of order *q*
#' (also called the “effective number of species”):
#'
#' \deqn{^{q}D \;=\;
#'    \left( \sum_{i} p_i^{\,q} \right)^{1/(1-q)}, \quad q \neq 1}
#'
#' For *q = 1* the formula is undefined; the limit is
#' \deqn{^{1}D = e^{H'}}.
#'
#' @param q Numeric order of diversity.  Common values:
#'   *0* (richness), *1* (exp(*H′*)), *2* (inverse Simpson).
#'
#' @return A **closure**: `hill_q(q)` returns a function that takes a
#'   vector of counts and yields the corresponding ^qD.  The returned
#'   function is vectorised over its input.
#'
#' @section References:
#' Hill, M. O. (1973) *Diversity and Evenness: A Unifying Notation and its
#' Consequences.* Ecology **54** (2), 427–432.
#'
#' @examples
#' hill1 <- hill_q(1)   # q = 1
#' hill1(c(5, 1, 1, 1))
#'
#' hill2 <- hill_q(2)   # q = 2, inverse-Simpson
#' hill2(c(5, 1, 1, 1))
#' @export
hill_q <- function(q) {
  force(q)
  function(cnt) {
    cnt <- cnt[cnt > 0]
    p   <- cnt / sum(cnt)
    if (q == 1) {
      exp(-sum(p * log(p)))
    } else {
      (sum(p ^ q))^(1 / (1 - q))
    }
  }
}



