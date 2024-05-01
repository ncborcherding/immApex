# test script for generate.sequences.R - testcases are NOT comprehensive!

test_that("generate.sequences works", {
  
  set.seed(42)
  t1.sequences <- generate.sequences(prefix.motif = "CAS",
                                     suffix.motif = "YF",
                                     number.of.sequences = 1000,
                                     min.length = 8,
                                     max.length = 16)
  
  expect_equal(
    t1.sequences,
    getdata("generate.sequences", "generate.sequences_T1")
  )
  
  set.seed(42)
  t2.sequences <- generate.sequences(number.of.sequences = 1000,
                                     min.length = 2,
                                     max.length = 16, 
                                     sequence.dictionary = c("A", "C", "T", "G"))
  
  expect_equal(
    t2.sequences,
    getdata("generate.sequences", "generate.sequences_T1")
  )
})