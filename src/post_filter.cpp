// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <tbb/concurrent_vector.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
using namespace Rcpp;
using namespace RcppParallel;

// Declaration for the edit_distance_threshold function (defined in edit_distance.cpp).
int edit_distance_threshold(std::string a, std::string b, int threshold);

// RcppParallel worker for post-filtering candidate pairs.
struct PostFilterWorker : public RcppParallel::Worker {
  const Rcpp::IntegerMatrix candidatePairs;
  const std::vector<std::string>& sequences;
  const std::vector<std::string>& vGenes;
  const std::vector<std::string>& jGenes;
  int threshold;
  bool filterV;
  bool filterJ;
  
  tbb::concurrent_vector< std::pair<int,int> > validPairs;
  
  PostFilterWorker(const Rcpp::IntegerMatrix& cand,
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
        validPairs.push_back(std::make_pair(idx1, idx2));
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::IntegerMatrix post_filter_candidates(const Rcpp::IntegerMatrix& candidatePairs,
                                             const std::vector<std::string>& sequences,
                                             const std::vector<std::string>& vGenes,
                                             const std::vector<std::string>& jGenes,
                                             int threshold,
                                             bool filterV,
                                             bool filterJ) {
  PostFilterWorker worker(candidatePairs, sequences, vGenes, jGenes, threshold, filterV, filterJ);
  parallelFor(0, candidatePairs.nrow(), worker);
  
  int nPairs = worker.validPairs.size();
  Rcpp::IntegerMatrix result(nPairs, 2);
  for (int i = 0; i < nPairs; i++) {
    result(i, 0) = worker.validPairs[i].first;
    result(i, 1) = worker.validPairs[i].second;
  }
  return result;
}
