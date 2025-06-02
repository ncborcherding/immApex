# test script for sequenceEncoder .R - testcases are NOT comprehensive!

# -------------------------------------------------------------------------
# 1.  ONE-HOT 
# -------------------------------------------------------------------------
test_that("one-hot encoding returns correct shape and content", {
  
  seqs <- c("ACD",   # len 3
            "W",     # len 1  -> will be padded
            "MNPQR") # len 5  -> will be truncated to max.length = 4
  
  res <- sequenceEncoder(seqs,
                         mode          = "onehot",
                         max.length    = 4,          # force trunc/pad
                         padding.symbol = ".",
                         verbose       = FALSE)
  
  expect_type(res,  "list")
  expect_named(res, c("cube", "flattened", "alphabet", "pad_token", "threads"),
               ignore.order = TRUE)
  
  # depth should be alphabet + padding = 21
  expect_equal(dim(res$cube), c(21, 4, length(seqs)))
  expect_equal(dim(res$flattened), c(length(seqs), 21 * 4))
  
  ## every position has *exactly one* 1 in the depth dimension
  one.count <- apply(res$cube, 3, sum)
  expect_equal(one.count, rep(4, length(seqs)))
  
  ## wrapper produces identical result
  res.w <- onehotEncoder(seqs, max.length = 4, verbose = FALSE)
  expect_equal(res$flattened, res.w$flattened)
})

# -------------------------------------------------------------------------
# 2.  PROPERTY (custom matrix) 
# -------------------------------------------------------------------------
test_that("property encoding with explicit matrix works and returns summary", {
  
  # toy property matrix: 20 AAs × 2 scales
  set.seed(42)
  prop.mat <- matrix(runif(40), nrow = 20,
                     dimnames = list(amino.acids, c("scale1", "scale2")))
  
  seqs <- c("AAAA", "CC", "EFG")          # different lengths
  
  res <- sequenceEncoder(seqs,
                         mode             = "property",
                         property.matrix  = prop.mat,
                         verbose          = FALSE)
  
  # depth = ncol(prop.mat) = 2
  expect_equal(dim(res$cube),       c(2, max(nchar(seqs)), length(seqs)))
  expect_equal(dim(res$flattened),  c(length(seqs), 2 * max(nchar(seqs))))
  expect_equal(dim(res$summary),    c(length(seqs), 2))
  
  ## manual means should match C++ summary
  manual.mean <- vapply(seqs, function(s) {
    aa <- strsplit(s, "")[[1]]
    colMeans(prop.mat[aa, , drop = FALSE])
  }, numeric(ncol(prop.mat)))
  
  expect_equal(manual.mean, res$summary, tolerance = 1e-12)
  
  ## wrapper equivalence
  res.w <- propertyEncoder(seqs,
                           property.matrix = prop.mat,
                           summary.fun     = "mean",
                           verbose         = FALSE)
  expect_equal(res$flattened, res.w$flattened)
})

# -------------------------------------------------------------------------
# 3.  PROPERTY via property.set  (skip if Peptides not installed) 
# -------------------------------------------------------------------------
test_that("property.set works when Peptides is installed", {
  
  skip_if_not_installed("Peptides")
  
  seqs <- c("ACDE", "W")
  
  res <- sequenceEncoder(seqs,
                         mode          = "property",
                         property.set  = c("hydrophobicity", "charge"),
                         summary.fun   = "",
                         verbose       = FALSE)
  
  expect_true(is.list(res))
  expect_null(res$summary)                   # no summary.fun requested
  expect_equal(dim(res$cube)[1],             2)  # two property scales
})

# -------------------------------------------------------------------------
# 4.  ERROR handling 
# -------------------------------------------------------------------------
test_that("informative errors are thrown for invalid input", {
  
  ## property mode without matrix or set
  expect_error(
    sequenceEncoder("AC", mode = "property", verbose = FALSE),
    "property.set` or `property.matrix` must be supplied"
  )
  
  ## non-existent property name
  if (requireNamespace("Peptides", quietly = TRUE)) {
    expect_error(
      sequenceEncoder("ACD",
                      mode = "property",
                      property.set = "not_a_real_property",
                      verbose = FALSE),
      "not_a_real_property"
    )
  }
  
  ## property.matrix wrong dimensions
  bad.mat <- matrix(1, nrow = 5, ncol = 1)
  expect_error(
    sequenceEncoder("AC",
                    mode = "property",
                    property.matrix = bad.mat,
                    verbose = FALSE),
    "rows"
  )
  
  ## zero sequences
  expect_error(sequenceEncoder(character(0), verbose = FALSE),
               "`sequences` is empty", fixed = TRUE)
})

# -------------------------------------------------------------------------
# 5.  Meta information 
# -------------------------------------------------------------------------
test_that("meta slots are coherent", {
  
  seqs <- c("AC", "DFG")
  res  <- sequenceEncoder(seqs, verbose = FALSE)
  
  expect_true(res$pad_token == ".")
  expect_true(res$alphabet[ res$pad_token ] == res$pad_token ||
                res$pad_token %in% res$alphabet)
  
  expect_true(res$threads >= 1)
})
