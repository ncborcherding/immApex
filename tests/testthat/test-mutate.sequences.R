# test script for mutate.sequences.R - testcases are NOT comprehensive!

test_that("mutate.sequences works", {
  
  sequences <- getdata("generate.sequences", "generate.sequences_T1")

  set.seed(42)
  mutate.default <- mutate.sequences(sequences)
  
  expect_equal(
    mutate.default,
    getdata("mutate.sequences", "mutate.sequences_default")
  )
  
  set.seed(42)
  mutate.mutation.rate <- mutate.sequences(sequences,
                                           mutation.rate = 0.2)
  
  expect_equal(
    mutate.mutation.rate,
    getdata("mutate.sequences", "mutate.sequences_mutation.rate")
  )
  
  set.seed(42)
  mutate.set.position <- mutate.sequences(sequences,
                                          position.start = 4,
                                          position.end = 8)
  
  expect_equal(
    mutate.set.position,
    getdata("mutate.sequences", "mutate.sequences_set.position")
  )
  
  set.seed(42)
  mutate.specific.mutation <- mutate.sequences(sequences,
                                               position.start = 4,
                                               position.end = 8,
                                               sequence.dictionary = c("V", "L", "I"))
  
  expect_equal(
    mutate.specific.mutation,
    getdata("mutate.sequences", "mutate.sequences_specific.mutation")
  )
  
  nt.sequences <- getdata("generate.sequences", "generate.sequences_T2")
  
  set.seed(42)
  mutate.nt <- mutate.sequences(nt.sequences,
                                sequence.dictionary = c("A", "C", "T", "G"))
  
  expect_equal(
    mutate.nt,
    getdata("mutate.sequences", "mutate.sequences_nucleotides")
  )
})

