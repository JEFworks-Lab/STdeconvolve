// Helpers for genes x pixels count matrices in compressed sparse column form
// (dgCMatrix). They replace several R-level passes over all the stored
// entries, each of which allocates vectors as long as the number of entries.

#include <Rcpp.h>
#include <algorithm>
#include <utility>
#include <vector>

// Per-row (gene) statistics over the columns (pixels) in `cols` (all columns
// if `cols` is empty):
//   sum       sum of the entries
//   detected  number of entries > 0
//   ssdev     sum of squared deviations from the row mean, over all columns
//             (only if `deviations`)
// [[Rcpp::export]]
Rcpp::List csc_row_stats(const Rcpp::IntegerVector& p, const Rcpp::IntegerVector& i,
                         const Rcpp::NumericVector& x, int nrow,
                         const Rcpp::LogicalVector& cols, bool deviations) {
  const int ncol = p.size() - 1;
  const bool all = cols.size() == 0;
  if (!all && cols.size() != ncol) Rcpp::stop("`cols` must have one value per column");
  Rcpp::NumericVector sum(nrow), detected(nrow), ssdev;
  std::vector<double> stored(nrow, 0.0);
  double n = 0.0;
  for (int j = 0; j < ncol; ++j) {
    if (!all && cols[j] != TRUE) continue;
    n += 1.0;
    for (int k = p[j]; k < p[j + 1]; ++k) {
      const int r = i[k];
      sum[r] += x[k];
      stored[r] += 1.0;
      if (x[k] > 0) detected[r] += 1.0;
    }
  }
  if (deviations) {
    ssdev = Rcpp::NumericVector(nrow);
    std::vector<double> mean(nrow);
    for (int r = 0; r < nrow; ++r) mean[r] = sum[r] / n;
    for (int j = 0; j < ncol; ++j) {
      if (!all && cols[j] != TRUE) continue;
      for (int k = p[j]; k < p[j + 1]; ++k) {
        const double d = x[k] - mean[i[k]];
        ssdev[i[k]] += d * d;
      }
    }
    // the entries that are not stored are zero
    for (int r = 0; r < nrow; ++r) ssdev[r] += (n - stored[r]) * mean[r] * mean[r];
  }
  return Rcpp::List::create(Rcpp::Named("sum") = sum,
                            Rcpp::Named("detected") = detected,
                            Rcpp::Named("ssdev") = ssdev,
                            Rcpp::Named("n") = n);
}

// Subset rows and columns in one pass. `row_map[r]` is the (1-based) row of
// row r in the result, or 0 to drop the row; `cols` are the columns to keep
// (all if empty), in their original order.
// [[Rcpp::export]]
Rcpp::List csc_subset(const Rcpp::IntegerVector& p, const Rcpp::IntegerVector& i,
                      const Rcpp::NumericVector& x,
                      const Rcpp::IntegerVector& row_map,
                      const Rcpp::LogicalVector& cols) {
  const int ncol = p.size() - 1;
  const bool all = cols.size() == 0;
  if (!all && cols.size() != ncol) Rcpp::stop("`cols` must have one value per column");
  bool ordered = true; // do the kept rows stay in their order?
  for (int r = 0, last = 0; r < row_map.size(); ++r) {
    if (row_map[r] == 0) continue;
    if (row_map[r] < last) { ordered = false; break; }
    last = row_map[r];
  }
  R_xlen_t nnz = 0;
  int nc = 0;
  for (int j = 0; j < ncol; ++j) {
    if (!all && cols[j] != TRUE) continue;
    ++nc;
    for (int k = p[j]; k < p[j + 1]; ++k) nnz += row_map[i[k]] != 0;
  }
  Rcpp::IntegerVector np(nc + 1), ni(nnz);
  Rcpp::NumericVector nx(nnz);
  std::vector<std::pair<int, double> > buf;
  R_xlen_t at = 0;
  int c = 0;
  for (int j = 0; j < ncol; ++j) {
    if (!all && cols[j] != TRUE) continue;
    const R_xlen_t begin = at;
    for (int k = p[j]; k < p[j + 1]; ++k) {
      const int r = row_map[i[k]];
      if (r == 0) continue;
      ni[at] = r - 1;
      nx[at] = x[k];
      ++at;
    }
    if (!ordered && at - begin > 1) {
      buf.resize(at - begin);
      for (R_xlen_t k = begin; k < at; ++k) buf[k - begin] = std::make_pair(ni[k], nx[k]);
      std::sort(buf.begin(), buf.end());
      for (R_xlen_t k = begin; k < at; ++k) {
        ni[k] = buf[k - begin].first;
        nx[k] = buf[k - begin].second;
      }
    }
    np[++c] = static_cast<int>(at);
  }
  return Rcpp::List::create(Rcpp::Named("p") = np, Rcpp::Named("i") = ni, Rcpp::Named("x") = nx);
}
