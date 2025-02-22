// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Banded dynamic programming edit distance with threshold.
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
