# test script for positional.encoder.R - testcases are NOT comprehensive!

test_that("positional.encode works", {
  
  set.seed(42)
  positional.encoding <- positional.encoder(number.of.sequences = 1000, 
                                            latent.dims = 64)
  
  expect_equal(
    positional.encoding,
    getdata("positional.encoding", "positional.encoding_values")
  )
  
})

