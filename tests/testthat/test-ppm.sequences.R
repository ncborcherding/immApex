# test script for PPM.sequences.R - testcases are NOT comprehensive!

test_that("PPM.sequences works", {
  
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
  
  set.seed(42)
  back.freq <- sample(1:1000, 20)
  back.freq <- back.freq/sum(back.freq)
  pwm.bf <- PPM.sequences(sequences,
                          max.length = 20,
                          convert.PWM = TRUE,
                          background.frequencies = back.freq)
  
  expect_equal(
    pwm.bf,
    getdata("ppm.sequences", "ppm.sequences_pwm.bf")
  )
  
  
})
