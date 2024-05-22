# test script for adjacency.matrix.R - testcases are NOT comprehensive!

test_that("adjacency.matrix works", {
  
  sequences <- getdata("generate.sequences", "generate.sequences_T1")

  adjacency.default <- adjacency.matrix(sequences)
  
  expect_equal(
    adjacency.default,
    getdata("adjacency.matrix", "adjacency.matrix_default")
  )
  
  adjacency.raw <- adjacency.matrix(sequences,
                                    normalize = FALSE)
  
  expect_equal(
    adjacency.raw,
    getdata("adjacency.matrix", "adjacency.matrix_unnormalized")
  )
  
  #TODO Add non standard sequence.dictionary
  
  
})
