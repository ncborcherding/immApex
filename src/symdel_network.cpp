// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <tbb/concurrent_vector.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

// ---------------------------
// Helper functions for deletion variants
// ---------------------------

// Recursively generate all combinations of indices to delete.
// "comb" holds the current combination (indices to delete)
// "variants" accumulates the deletion variants of the string s.
void generate_combinations(const std::string &s, int k, int start,
                           std::vector<int>& comb, std::vector<std::string>& variants) {
  if (k == 0) {
    std::string variant;
    int n = s.size();
    int combIndex = 0;
    for (int i = 0; i < n; i++) {
      if (combIndex < comb.size() && comb[combIndex] == i) {
        combIndex++;
      } else {
        variant.push_back(s[i]);
      }
    }
    variants.push_back(variant);
    return;
  }
  for (int i = start; i <= s.size() - k; i++) {
    comb.push_back(i);
    generate_combinations(s, k - 1, i + 1, comb, variants);
    comb.pop_back();
  }
}

// Generate all deletion variants for string s with up to d deletions.
std::vector<std::string> generate_deletion_variants(const std::string &s, int d) {
  std::vector<std::string> variants;
  for (int k = 0; k <= d; k++) {
    std::vector<int> comb;
    generate_combinations(s, k, 0, comb, variants);
  }
  return variants;
}

// ---------------------------
// Symmetric Deletion Lookup
// ---------------------------

// [[Rcpp::export]]
IntegerMatrix symmetric_deletion_lookup_cpp(std::vector<std::string> sequences, int threshold) {
  int n = sequences.size();
  std::vector<int> lengths(n);
  for (int i = 0; i < n; i++) {
    lengths[i] = sequences[i].size();
  }
  
  // Build a hash map: deletion variant string -> vector of sequence indices.
  std::unordered_map<std::string, std::vector<int>> variant_map;
  for (int i = 0; i < n; i++) {
    std::vector<std::string> variants = generate_deletion_variants(sequences[i], threshold);
    for (const std::string &var : variants) {
      variant_map[var].push_back(i);
    }
  }
  
  // Use a set to collect candidate pairs without duplicates.
  std::set< std::pair<int,int> > candidatePairs;
  for (int i = 0; i < n; i++) {
    std::vector<std::string> variants = generate_deletion_variants(sequences[i], threshold);
    for (const std::string &var : variants) {
      if (variant_map.find(var) != variant_map.end()) {
        const std::vector<int> &indices = variant_map[var];
        for (int j : indices) {
          if (i == j) continue;
          // Basic length prefilter.
          if (std::abs(lengths[i] - lengths[j]) > threshold) continue;
          int a = std::min(i, j);
          int b = std::max(i, j);
          candidatePairs.insert(std::make_pair(a, b));
        }
      }
    }
  }
  
  // Convert candidate pairs to an R IntegerMatrix (1-indexed).
  int numPairs = candidatePairs.size();
  IntegerMatrix out(numPairs, 2);
  int row = 0;
  for (auto it = candidatePairs.begin(); it != candidatePairs.end(); ++it) {
    out(row, 0) = it->first + 1;
    out(row, 1) = it->second + 1;
    row++;
  }
  return out;
}

// ---------------------------
// Banded Edit Distance with Threshold
// ---------------------------

// Compute edit distance with a given threshold using banded dynamic programming.
// Returns a value > threshold if the edit distance exceeds the threshold.
int edit_distance_threshold_cpp(const std::string &a, const std::string &b, int threshold) {
  int n = a.size();
  int m = b.size();
  if (std::abs(n - m) > threshold) return threshold + 1;
  
  std::vector<int> prev(m+1), curr(m+1);
  for (int j = 0; j <= m; j++) {
    prev[j] = j;
  }
  
  for (int i = 1; i <= n; i++) {
    curr[0] = i;
    int j_start = std::max(1, i - threshold);
    int j_end = std::min(m, i + threshold);
    for (int j = 1; j < j_start; j++) {
      curr[j] = threshold + 1;
    }
    for (int j = j_start; j <= j_end; j++) {
      int cost = (a[i-1] == b[j-1]) ? 0 : 1;
      curr[j] = std::min({ prev[j] + 1, curr[j-1] + 1, prev[j-1] + cost });
    }
    for (int j = j_end+1; j <= m; j++) {
      curr[j] = threshold + 1;
    }
    prev.swap(curr);
  }
  return prev[m];
}

// [[Rcpp::export]]
int edit_distance_threshold(std::string a, std::string b, int threshold) {
  return edit_distance_threshold_cpp(a, b, threshold);
}

// ---------------------------
// RcppParallel Worker for Post-filtering Candidate Pairs
// (This worker computes the edit distance and retains it as the edge weight.)
// ---------------------------

// Define a structure to hold a triple: from, to, and weight.
struct Triple {
  int from;
  int to;
  int weight;
};

struct PostFilterWorker : public Worker {
  const IntegerMatrix candidatePairs;
  const std::vector<std::string>& sequences;
  const std::vector<std::string>& vGenes;
  const std::vector<std::string>& jGenes;
  int threshold;
  bool filterV;
  bool filterJ;
  
  tbb::concurrent_vector<Triple> validTriples;
  
  PostFilterWorker(const IntegerMatrix& cand,
                   const std::vector<std::string>& seq,
                   const std::vector<std::string>& vG,
                   const std::vector<std::string>& jG,
                   int thresh, bool filtV, bool filtJ)
    : candidatePairs(cand), sequences(seq), vGenes(vG), jGenes(jG),
      threshold(thresh), filterV(filtV), filterJ(filtJ) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      int idx1 = candidatePairs(i, 0);
      int idx2 = candidatePairs(i, 1);
      int index1 = idx1 - 1;
      int index2 = idx2 - 1;
      
      if (filterV && (vGenes[index1] != vGenes[index2])) continue;
      if (filterJ && (jGenes[index1] != jGenes[index2])) continue;
      
      int d = edit_distance_threshold(sequences[index1], sequences[index2], threshold);
      if (d <= threshold) {
        validTriples.push_back({idx1, idx2, d});
      }
    }
  }
};

// [[Rcpp::export]]
DataFrame post_filter_candidates(const IntegerMatrix& candidatePairs,
                                 const std::vector<std::string>& sequences,
                                 const std::vector<std::string>& vGenes,
                                 const std::vector<std::string>& jGenes,
                                 int threshold,
                                 bool filterV,
                                 bool filterJ) {
  PostFilterWorker worker(candidatePairs, sequences, vGenes, jGenes, threshold, filterV, filterJ);
  parallelFor(0, candidatePairs.nrow(), worker);
  
  int nTriples = worker.validTriples.size();
  std::vector<int> from(nTriples), to(nTriples), weight(nTriples);
  for (int i = 0; i < nTriples; i++) {
    from[i] = worker.validTriples[i].from;
    to[i] = worker.validTriples[i].to;
    weight[i] = worker.validTriples[i].weight;
  }
  return DataFrame::create(Named("from") = from,
                           Named("to") = to,
                           Named("weight") = weight);
}
