#include <Rcpp.h>
using namespace Rcpp;

/* ------------------------------------------------------------------------- */
/*  Levenshtein distance with early exit when                */
/*  current row’s minimum already exceeds user threshold.     */
/* ------------------------------------------------------------------------- */
static inline int lv_threshold(const std::string &s,
                               const std::string &t,
                               int thresh)
{
  const int n = s.size(), m = t.size();
  if (std::abs(n - m) > thresh)        // impossible to be ≤ thresh
    return thresh + 1;
  
  std::vector<int> prev(m + 1), cur(m + 1);
  std::iota(prev.begin(), prev.end(), 0);
  
  for (int i = 1; i <= n; ++i) {
    cur[0]   = i;
    int row_min = cur[0];
    
    for (int j = 1; j <= m; ++j) {
      const int cost = (s[i - 1] == t[j - 1]) ? 0 : 1;
      cur[j] = std::min({ prev[j] + 1,
                        cur[j - 1] + 1,
                        prev[j - 1] + cost });
      row_min = std::min(row_min, cur[j]);
    }
    if (row_min > thresh) return thresh + 1;
    std::swap(prev, cur);
  }
  return prev[m];
}

/* ------------------------------------------------------------------------- */
/*  Vector-wise fast edit-distance → edge list                                */
/* ------------------------------------------------------------------------- */
// [[Rcpp::export]]
DataFrame fast_edge_list(CharacterVector           seqs,
                         int                       thresh,
                         Nullable<CharacterVector> ids = R_NilValue)
{
  const int n = seqs.size();
  std::vector<std::string> s(n);
  for (int i = 0; i < n; ++i)
    s[i] = as<std::string>(seqs[i]);
  
  /* ---------------- handle labels ---------------- */
  std::vector<std::string> lbl;
  
  if (ids.isNull()) {                       // nothing supplied
    lbl.reserve(n);
    for (int i = 0; i < n; ++i)
      lbl.emplace_back(std::to_string(i + 1));
  } else {                                  // convert the vector
    CharacterVector ids_cv(ids);            // extract from Nullable
    lbl = as< std::vector<std::string> >(ids_cv);
    if ((int)lbl.size() != n)
      stop("Length of 'ids' must equal length of 'seqs'.");
  }
  
  /* --------------- rest unchanged --------------- */
  std::vector<std::string> out_from, out_to;
  std::vector<int>         out_dist;
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<std::string> local_from, local_to;
  std::vector<int>         local_dist;
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 32)
#endif
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (std::abs((int)s[i].size() - (int)s[j].size()) > thresh) continue;
      int d = lv_threshold(s[i], s[j], thresh);
      if (d > thresh) continue;
      
      local_from.emplace_back(lbl[i]);
      local_to  .emplace_back(lbl[j]);
      local_dist.emplace_back(d);
    }
  }
  
#ifdef _OPENMP
#pragma omp critical
#endif
{
  out_from.insert(out_from.end(), local_from.begin(), local_from.end());
  out_to  .insert(out_to  .end(), local_to  .begin(), local_to  .end());
  out_dist.insert(out_dist.end(),local_dist.begin(),local_dist.end());
}
}

return DataFrame::create(_["from"] = out_from,
                         _["to"]   = out_to,
                         _["dist"] = out_dist,
                         _["stringsAsFactors"] = false);
}