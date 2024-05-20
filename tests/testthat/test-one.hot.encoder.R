# test script for one.hot.encoder.R - testcases are NOT comprehensive!

test_that("one.hot.encoder works", {
  
  sequences <- getdata("generate.sequences", "generate.sequences_T1")

  ppm.default <- PPM.sequences(sequences)
  
  expect_equal(
    ppm.default,
    getdata("ppm.sequences", "ppm.sequences_default")
  )
  
  pwm.default <- PPM.sequences(sequences,
                               convert.PWM = TRUE)
  
  expect_equal(
    pwm.default,
    getdata("ppm.sequences", "ppm.sequences_pwm.default")
  )
  
  
})
