# test script for geometric.encoder.R - testcases are NOT comprehensive!

test_that("geometric.encoder works", {
  
  sequences <- getdata("generate.sequences", "generate.sequences_T1")

  geom.default <- geometric.encoder(sequences,
                                    theta = pi)
  
  expect_equal(
    geom.default,
    getdata("geometric.encoder", "geometric.encoder_matrix")
  )
  
  geom.BLOSUM45 <- geometric.encoder(sequences, 
                                     method.to.use = "BLOSUM45",
                                     theta = pi/3)
  expect_equal(
    geom.BLOSUM45,
    getdata("geometric.encoder", "geometric.encoder_BLOSUM45matrix")
  )
  
  geom.PAM30 <- geometric.encoder(sequences, 
                                  method.to.use = "PAM30",
                                  theta = pi/2)
  expect_equal(
    geom.PAM30,
    getdata("geometric.encoder", "geometric.encoder_PAM30matrix")
  )
  
  
})
