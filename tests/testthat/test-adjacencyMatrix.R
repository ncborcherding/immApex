# test script for adjacencyMatrix.R - testcases are NOT comprehensive!

test_that("adjacencyMatrix works", {
  
  sequences <- getdata("generateSequences", "generateSequences_T1")

  adjacency.default <- adjacencyMatrix(sequences)
  
  expect_equal(
    adjacency.default,
    getdata("adjacencyMatrix", "adjacencyMatrix_default")
  )
  
  adjacency.raw <- adjacencyMatrix(sequences,
                                    normalize = FALSE)
  
  expect_equal(
    adjacency.raw,
    getdata("adjacencyMatrix", "adjacencyMatrix_unnormalized")
  )
  
  #TODO Add non standard sequence.dictionary
  
  
})
