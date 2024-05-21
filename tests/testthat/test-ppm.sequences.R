# test script for probability.matrix.R - testcases are NOT comprehensive!

test_that("probability.matrix works", {
  
  sequences <- getdata("generate.sequences", "generate.sequences_T1")

  ppm.default <- probability.matrix(sequences)
  
  expect_equal(
    ppm.default,
    getdata("ppm.sequences", "ppm.sequences_default")
  )
  
  pwm.default <- probability.matrix(sequences,
                                    convert.PWM = TRUE)
  
  expect_equal(
    pwm.default,
    getdata("ppm.sequences", "ppm.sequences_pwm.default")
  )
  
  set.seed(42)
  back.freq <- sample(1:1000, 20)
  back.freq <- back.freq/sum(back.freq)
  pwm.bf <- probability.matrix(sequences,
                               max.length = 20,
                               convert.PWM = TRUE,
                               background.frequencies = back.freq)
  
  expect_equal(
    pwm.bf,
    getdata("ppm.sequences", "ppm.sequences_pwm.bf")
  )
  
  
})
