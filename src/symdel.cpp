// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// Helper: recursively generate all combinations of indices to delete.
// "comb" holds the current combination (indices to delete)
// "variants" accumulates the deletion variants of the string s.
void generate_combinations(const std::string &s, int k, int start, std::vector<int>& comb, std::vector<std::string>& variants) {
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
