# test script for geometricEncoder.R - testcases are NOT comprehensive!

mock_identity_matrix <- diag(20)
rownames(mock_identity_matrix) <- amino.acids
colnames(mock_identity_matrix) <- amino.acids # Good practice

mock_simple_matrix <- matrix(1:400, 20, 20,
                             dimnames = list(amino.acids, NULL))

immapex_blosum.pam.matrices <- list(
  BLOSUM62 = mock_identity_matrix,
  PAM250 = mock_simple_matrix
)

test_that("Basic encoding with identity matrix and no rotation is correct", {
  
  emb_A <- geometricEncoder(c("A"), theta = 2 * pi)
  emb_AR <- geometricEncoder(c("AR"), theta = 2 * pi)
  
  expected_A <- c(1, rep(0, 19))
  expected_AR <- c(0.5, 0.5, rep(0, 18))
  
  expect_equal(as.vector(emb_A), expected_A)
  expect_equal(as.vector(emb_AR), expected_AR)
  expect_equal(dim(emb_AR), c(1, 20))
})

test_that("Encoding with 90-degree rotation (pi/2) is correct", {
  
  emb_A <- geometricEncoder("A", theta = pi / 2)
  emb_R <- geometricEncoder("R", theta = pi / 2)
  
  expected_A <- c(0, -1, rep(0, 18))
  expected_R <- c(1, 0, rep(0, 18))
  
  expect_equal(as.vector(emb_A), expected_A, tolerance = 1e-9)
  expect_equal(as.vector(emb_R), expected_R, tolerance = 1e-9)
})

test_that("Matrix selection logic works correctly", {
  emb_default <- geometricEncoder("A", theta = 0)
  expect_equal(as.vector(emb_default), as.vector(mock_identity_matrix[1, ]))
  
  # Fetches PAM250 from the list
  emb_pam <- geometricEncoder("A", method = "PAM250", theta = 0)
  expect_equal(as.vector(emb_pam), as.vector(mock_simple_matrix[1, ]))
  
  # New Feature: Correctly uses a matrix passed directly to `method`
  custom_mat <- matrix(rep(1:20, each = 20), 20, 20, dimnames = list(amino.acids, NULL))
  emb_custom <- geometricEncoder("A", method = custom_mat, theta = 0)
  expect_equal(as.vector(emb_custom), custom_mat[1, ])
})

test_that("Output dimensions are correct for multiple sequences", {
  seqs <- c("A", "R", "AR", "CAS")
  emb <- geometricEncoder(seqs)
  expect_equal(dim(emb), c(length(seqs), 20))
  expect_true(is.numeric(emb))
})

test_that("Error handling for invalid inputs is robust", {
  # Gracefully handles empty character vector input
  expect_equal(dim(geometricEncoder(character(0))), c(0, 20))
  
  # Rejects NA and empty strings with a specific error
  err_msg <- "NA or empty strings are not allowed"
  expect_error(geometricEncoder(c("A", NA)), err_msg)
  expect_error(geometricEncoder(c("A", "")), err_msg)
  
  # Rejects non-canonical amino acids with a specific error
  expect_error(geometricEncoder("AX"), "Non-canonical amino acid\\(s\\) found: X")
  expect_error(geometricEncoder("BZX"), "Non-canonical amino acid\\(s\\) found: B, Z, X")
  
  # Rejects invalid method key
  expect_error(geometricEncoder("A", method = "FOOBAR"), "Cannot find matrix for method 'FOOBAR'.")
  
  # Rejects invalid matrices passed directly to `method`
  bad_mat_dim <- matrix(1, 5, 5)
  bad_mat_rownames <- matrix(1, 20, 20)
  expect_error(geometricEncoder("A", method = bad_mat_dim), "must be 20x20")
  expect_error(geometricEncoder("A", method = bad_mat_rownames), "with amino acid rownames")
})
