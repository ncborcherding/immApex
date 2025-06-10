# test script for variationalSequences.R - testcases are NOT comprehensive!

set.seed(42)
sequences <-  generateSequences(prefix.motif = "CAS",
                                suffix.motif = "YF",
                                number.of.sequences = 100,
                                min.length = 8,
                                max.length = 16)
                         
test_that("Test for valid input sequence length", {
        expect_error(variationalSequences(input.sequences = character(0)),
                     "input.sequences must have at least one sequence.")
})
      
test_that("Test for valid mode", {
        expect_error(variationalSequences(input.sequences = sequences,
                                          mode = "invalidEncoder", 
                                          verbose = FALSE),
                    regexp = "onehot")
})
      
test_that("Test for valid optimizer", {
        expect_error(variationalSequences(input.sequences = sequences,
                                          optimizer = "invalidOptimizer"),
                     "Please select a compatible optimizer function in the Keras R implementation.")
})
      
test_that("Test for correct output type", {
        result <- variationalSequences(input.sequences = sequences, 
                                       sequence.dictionary = amino.acids[1:20])
        expect_true(is.vector(result))
})
      
test_that("Test for correct sequence generation", {
        number_of_sequences <- 5
        result <- variationalSequences(input.sequences = sequences, 
                                       number.of.sequences = number_of_sequences, 
                                       verbose = FALSE)
        expect_equal(length(result), number_of_sequences)
})

