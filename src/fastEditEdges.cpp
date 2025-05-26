// ── src/fastEditEdges_group.cpp ─────────────────────────────────────────────
#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;

/*---------------------------- helper from before ---------------------------*/
static inline int lv_threshold(const std::string &s,
                               const std::string &t,
                               int thresh)
{
  const int n = s.size(), m = t.size();
  if (std::abs(n - m) > thresh) return thresh + 1;
  
  std::vector<int> prev(m + 1), cur(m + 1);
  std::iota(prev.begin(), prev.end(), 0);
  
  for (int i = 1; i <= n; ++i) {
    cur[0] = i;
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
  return prev.back();
}

/*-------------------------------------------------------------------------*//**
 *  Fast edit-distance edge list with optional V / J matching.
 *
 *  @param seqs     Character vector of sequences.
 *  @param thresh   Integer  ≥1 (“absolute”) **or** 0–1 (“relative”).
 *  @param v_gene   (optional) Character vector; recycled if length 1.
 *  @param j_gene   (optional) Character vector; recycled if length 1.
 *  @param match_v  If TRUE keep only pairs with identical V gene.
 *  @param match_j  If TRUE keep only pairs with identical J gene.
 *  @param ids      (optional) labels; default 1:n
 *
 *  @return data.frame 〈from,to,dist〉
 */
// [[Rcpp::export]]
DataFrame fast_edge_list(CharacterVector           seqs,
                         double                    thresh          = 2.0,
                         Nullable<CharacterVector> v_gene = R_NilValue,
                         Nullable<CharacterVector> j_gene = R_NilValue,
                         bool                      match_v         = false,
                         bool                      match_j         = false,
                        Nullable<CharacterVector> ids    = R_NilValue)
{
  const int n = seqs.size();
  if (n < 2) stop("Need at least two sequences.");
  
  /* ---------- convert input vectors to std::string for speed ------------ */
  std::vector<std::string> s(n);
  for (int i = 0; i < n; ++i) s[i] = as<std::string>(seqs[i]);
  
  /* ---------- V / J annotation (optional) ------------------------------- */
  std::vector<std::string> v(n), j(n);
  if (match_v) {
    if (v_gene.isNull())
      stop("match_v = TRUE but v_gene is missing");
    CharacterVector vv(v_gene);
    if (vv.size() == 1) vv = CharacterVector(n, vv[0]);   // recycle scalar
    if (vv.size() != n) stop("v_gene length must equal length(seqs)");
    v = as< std::vector<std::string> >(vv);
  }
  if (match_j) {
    if (j_gene.isNull())
      stop("match_j = TRUE but j_gene is missing");
    CharacterVector jj(j_gene);
    if (jj.size() == 1) jj = CharacterVector(n, jj[0]);
    if (jj.size() != n) stop("j_gene length must equal length(seqs)");
    j = as< std::vector<std::string> >(jj);
  }
  
  /* ---------- labels ---------------------------------------------------- */
  std::vector<std::string> lbl;
  if (ids.isNull()) {
    lbl.reserve(n);
    for (int i = 0; i < n; ++i) lbl.emplace_back(std::to_string(i + 1));
  } else {
    CharacterVector ids_cv(ids);
    if (ids_cv.size() != n)
      stop("ids length must equal length(seqs)");
    lbl = as< std::vector<std::string> >(ids_cv);
  }
  
  /* ---------- output containers (thread-safe aggregation) --------------- */
  std::vector<std::string> out_from, out_to;
  std::vector<int>         out_dist;
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<std::string> local_from, local_to;
  std::vector<int>         local_dist;
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic,32)
#endif
  for (int i = 0; i < n; ++i) {
    const int len_i = s[i].size();
    
    for (int j2 = i + 1; j2 < n; ++j2) {
      
      /* ---- V/J filters ---------------------------------------------- */
      if (match_v && v[i] != v[j2]) continue;
      if (match_j && j[i] != j[j2]) continue;
      
      /* ---- length-difference pre-filter ----------------------------- */
      const int len_j = s[j2].size();
      int max_dist =
        (thresh <= 1.0)
        ? static_cast<int>(std::ceil(thresh * std::max(len_i, len_j)))
        : static_cast<int>(std::round(thresh));
      
      if (std::abs(len_i - len_j) > max_dist) continue;
      
      /* ---- edit distance ------------------------------------------- */
      int d = lv_threshold(s[i], s[j2], max_dist);
      if (d > max_dist) continue;
      
      /* ---- store ---------------------------------------------------- */
      local_from.emplace_back(lbl[i]);
      local_to  .emplace_back(lbl[j2]);
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
}  // end parallel section

return DataFrame::create(_["from"] = out_from,
                         _["to"]   = out_to,
                         _["dist"] = out_dist,
                         _["stringsAsFactors"] = false);
}
