# test script for probabilityMatrix .R - testcases are NOT comprehensive!

test_that("probabilityMatrix  works", {
  
  sequences <- getdata("generateSequences", "generateSequences_T1")

  ppm.default <- probabilityMatrix (sequences)
  
  expect_equal(
    ppm.default,
    getdata("probabilityMatrix", "ppm.sequences_default")
  )
  
  pwm.default <- probabilityMatrix (sequences,
                                    convert.PWM = TRUE)
  
  expect_equal(
    pwm.default,
    getdata("probabilityMatrix", "ppm.sequences_pwm.default")
  )
  
  set.seed(42)
  back.freq <- sample(1:1000, 20)
  back.freq <- back.freq/sum(back.freq)
  pwm.bf <- probabilityMatrix (sequences,
                               max.length = 20,
                               convert.PWM = TRUE,
                               background.frequencies = back.freq)
  
  expect_equal(
    pwm.bf,
    getdata("probabilityMatrix", "ppm.sequences_pwm.bf")
  )
  
  #TODO Add non standard sequence.dictionary
})
