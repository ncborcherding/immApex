# test script for geometricEncoder.R - testcases are NOT comprehensive!

test_that("geometricEncoder works", {
  
  sequences <- getdata("generateSequences", "generateSequences_T1")

  geom.default <- geometricEncoder(sequences,
                                    theta = pi)
  
  expect_equal(
    geom.default,
    getdata("geometricEncoder", "geometricEncoder_matrix")
  )
  
  geom.BLOSUM45 <- geometricEncoder(sequences, 
                                     method = "BLOSUM45",
                                     theta = pi/3)
  expect_equal(
    geom.BLOSUM45,
    getdata("geometricEncoder", "geometricEncoder_BLOSUM45matrix")
  )
  
  geom.PAM30 <- geometricEncoder(sequences, 
                                  method = "PAM30",
                                  theta = pi/2)
  expect_equal(
    geom.PAM30,
    getdata("geometricEncoder", "geometricEncoder_PAM30matrix")
  )
  
  
})
