# test script for one.hot.encoder.R - testcases are NOT comprehensive!

test_that("one.hot.encoder works", {
  
  sequences <- getdata("generate.sequences", "generate.sequences_T1")

  ohe.default <- one.hot.encoder(sequences)
  
  expect_equal(
    ohe.default,
    getdata("ohe.encoder", "one.hot.encoder_matrix")
  )
  
  ohe.padded <- one.hot.encoder(sequences,
                                max.length = 40)
  
  expect_equal(
    ohe.padded,
    getdata("ohe.encoder", "one.hot.encoder_padded.matrix")
  )
  
  ohe.array <- one.hot.encoder(sequences,
                               convert.to.matrix = FALSE)
  
  expect_equal(
    ohe.array,
    getdata("ohe.encoder", "one.hot.encoder_array")
  )
  
  ohe.padded.array <- one.hot.encoder(sequences,
                               max.length = 40,
                               convert.to.matrix = FALSE)
  
  expect_equal(
    ohe.padded.array,
    getdata("ohe.encoder", "one.hot.encoder_padded.array")
  )
  
  nt.sequences <- getdata("generate.sequences", "generate.sequences_T2")
  
  ohe.nt <- one.hot.encoder(nt.sequences,
                            sequence.dictionary = c("A", "C", "T", "G"))
  
  expect_equal(
    ohe.nt,
    getdata("ohe.encoder", "one.hot.encoder_nucleotide.matrix")
  )
  
})

#TODO Add motif testing