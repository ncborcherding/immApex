# test script for tokenize.sequences.R - testcases are NOT comprehensive!

test_that("tokenize.sequences works", {
  
  sequences <- getdata("generate.sequences", "generate.sequences_T1")

  token.default <- tokenize.sequences(sequences)
  
  expect_equal(
    token.default,
    getdata("tokenize.sequences", "tokenize.sequences_default")
  )
  
  token.noAdded <- tokenize.sequences(sequences,
                                      add.startstop = FALSE)
  
  expect_equal(
    token.noAdded,
    getdata("tokenize.sequences", "tokenize.sequences_noStartStop")
  )
  
  token.matrix <- tokenize.sequences(sequences,
                                     convert.to.matrix = TRUE)
  
  expect_equal(
    token.noAdded,
    getdata("tokenize.sequences", "tokenize.sequences_matrix")
  )
  
  
})

