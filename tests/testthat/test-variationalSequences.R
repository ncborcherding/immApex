# test script for variationalSequences.R - testcases are NOT comprehensive!

test_that("variationalSequences works", {
  
  sequences <- getdata("generateSequences", "generateSequences_T1")
  
  set.seed(42)
  new.sequences <- variationalSequences(sequences, 
                                        number.of.sequences = 10)
  
  expect_equal(
    length(new.sequences),
    10
  )
  expect_equal(
    substring(new.sequences[1], 1,3),
      "CAS"
  )
  
  new.sequences <- variationalSequences(sequences, 
                                        number.of.sequences = 10,
                                        encoder.function = "propertyEncoder",
                                        aa.method.to.use = "kideraFactors",
                                        call.threshold = 1)
  
  expect_equal(
    length(new.sequences),
    10
  )
  expect_equal(
    substring(new.sequences[1], 1,3),
    "CAS"
  )
  
})