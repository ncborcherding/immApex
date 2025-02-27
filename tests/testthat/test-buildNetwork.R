# test script for buildNetwork.R - testcases are NOT comprehensive!
library(igraph)

test_that("buildNetwork returns an igraph object and correct vertex count", {
  seqs <- c("AAA", "AAB", "ABA", "ABB")
  g <- buildNetwork(seqs, threshold = 1)
  expect_true(inherits(g, "igraph"))
  expect_equal(vcount(g), length(seqs))
})

test_that("Edge weights are computed correctly", {
  # For these sequences, "AAA" and "AAB" differ by 1, so we expect an edge with weight 1.
  seqs <- c("AAA", "AAB", "BBB")
  g <- buildNetwork(seqs, threshold = 1)
  # Check that at least one edge has weight 1.
  edge_weights <- E(g)$weight
  expect_true(any(edge_weights == 1))
  # Also, check that all weights are <= threshold 
  if (ecount(g) > 0) {
    expect_true(all(edge_weights <= 1))
  }
})

test_that("Function works with data frame input and retains gene annotations", {
  df <- data.frame(
    sequence = c("AAA", "AAB", "ABA", "ABB"),
    v = c("V1", "V1", "V2", "V1"),
    j = c("J1", "J1", "J2", "J1"),
    stringsAsFactors = FALSE
  )
  # With filtering enabled, only pairs with matching v.gene and j.gene should be connected.
  g <- buildNetwork(df, threshold = 1, filter.v = TRUE, filter.j = TRUE)
  expect_true(inherits(g, "igraph"))
  # Check that vertex attributes match the input.
  expect_equal(V(g)$sequence, df$sequence)
  expect_equal(V(g)$v.gene, df[["v.gene"]])
  expect_equal(V(g)$j.gene, df[["j.gene"]])
})

test_that("buildNetwork returns an empty edge set when no pairs qualify", {
  seqs <- c("AAA", "AAB", "ABA")
  g <- buildNetwork(seqs, threshold = 0)
  expect_equal(ecount(g), 0)
})
