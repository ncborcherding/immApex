# test script for positionalEncoder.R - testcases are NOT comprehensive!

test_that("positionalEncoder works", {
  
  set.seed(42)
  positional.encoding <- positionalEncoder(number.of.sequences = 1000, 
                                            latent.dims = 64)
  
  expect_equal(
    positional.encoding,
    getdata("positional.encoding", "positional.encoding_values")
  )
  
})

