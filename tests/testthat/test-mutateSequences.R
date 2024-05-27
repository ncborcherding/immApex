# test script for mutateSequences.R - testcases are NOT comprehensive!

test_that("mutateSequences works", {
  
  sequences <- getdata("generateSequences", "generateSequences_T1")

  set.seed(42)
  mutate.default <- mutateSequences(sequences)
  
  expect_equal(
    mutate.default,
    getdata("mutateSequences", "mutateSequences_default")
  )
  
  set.seed(42)
  mutate.mutation.rate <- mutateSequences(sequences,
                                           mutation.rate = 0.2)
  
  expect_equal(
    mutate.mutation.rate,
    getdata("mutateSequences", "mutateSequences_mutation.rate")
  )
  
  set.seed(42)
  mutate.set.position <- mutateSequences(sequences,
                                          position.start = 4,
                                          position.end = 8)
  
  expect_equal(
    mutate.set.position,
    getdata("mutateSequences", "mutateSequences_set.position")
  )
  
  set.seed(42)
  mutate.specific.mutation <- mutateSequences(sequences,
                                               position.start = 4,
                                               position.end = 8,
                                               sequence.dictionary = c("V", "L", "I"))
  
  expect_equal(
    mutate.specific.mutation,
    getdata("mutateSequences", "mutateSequences_specific.mutation")
  )
  
  nt.sequences <- getdata("generateSequences", "generateSequences_T2")
  
  set.seed(42)
  mutate.nt <- mutateSequences(nt.sequences,
                                sequence.dictionary = c("A", "C", "T", "G"))
  
  expect_equal(
    mutate.nt,
    getdata("mutateSequences", "mutateSequences_nucleotides")
  )
})

