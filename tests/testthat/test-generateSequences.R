# test script for generateSequences.R - testcases are NOT comprehensive!

test_that("generateSequences works", {
    
    set.seed(42)
    t1.sequences <- generateSequences(prefix.motif = "CAS",
                                      suffix.motif = "YF",
                                      number.of.sequences = 1000,
                                      min.length = 8,
                                      max.length = 16)
    
    expect_equal(
      t1.sequences,
      getdata("generateSequences", "new_generateSequences_T1")
    )
    
    set.seed(42)
    t2.sequences <- generateSequences(number.of.sequences = 1000,
                                      min.length = 2,
                                      max.length = 16, 
                                      sequence.dictionary = c("A", "C", "T", "G"))
    
    expect_equal(
      t2.sequences,
      getdata("generateSequences", "generateSequences_T2")
    )

})

