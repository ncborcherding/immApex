# test script for tokenizeSequences.R - testcases are NOT comprehensive!

test_that("tokenizeSequences works", {
  
  sequences <- getdata("generateSequences", "generateSequences_T1")

  token.default <- tokenizeSequences(sequences)
  
  expect_equal(
    token.default,
    getdata("tokenizeSequences", "tokenizeSequences_default")
  )
  
  token.noAdded <- tokenizeSequences(sequences,
                                      add.startstop = FALSE)
  
  expect_equal(
    token.noAdded,
    getdata("tokenizeSequences", "tokenizeSequences_noStartStop")
  )
  
  token.matrix <- tokenizeSequences(sequences,
                                     convert.to.matrix = TRUE)
  
  expect_equal(
    token.matrix,
    getdata("tokenizeSequences", "tokenizeSequences_matrix")
  )
  
  
})

