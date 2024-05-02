# test script for property.encoder.R - testcases are NOT comprehensive!

test_that("property.encoder works", {
  
  sequences <- getdata("generate.sequences", "generate.sequences_T1")
  
  #Return Matrix
  af.matrix <- property.encoder(sequences,
                                method.to.use = "atchleyFactors",
                                convert.to.matrix = TRUE)
  
  expect_equal(
    af.matrix,
    getdata("property.encoder", "property.encoder_AtchleyFactors_matrix")
  )
  
  #Return Array
  kf.array <- property.encoder(sequences,
                               method.to.use = "kideraFactors",
                               convert.to.matrix = FALSE)
  
  expect_equal(
    kf.array,
    getdata("property.encoder", "property.encoder_KideraFactors_array")
  )
  
  #Padded Matrix
  fasgai.matrix <- property.encoder(sequences,
                                    max.length = 40,
                                    method.to.use = "FASGAI",
                                    convert.to.matrix = FALSE)
  
  expect_equal(
    fasgai.matrix,
    getdata("property.encoder", "property.encoder_FASGAI_matrix")
  )
  
  #Padded Array
  vhse.array <- property.encoder(sequences,
                                 max.length = 40,
                                 method.to.use = "VHSE",
                                 convert.to.matrix = FALSE)
  
  expect_equal(
    vhse.array,
    getdata("property.encoder", "property.encoder_VSHE_array")
  )
  
})