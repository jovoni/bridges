// Prototype: G and B with length-saturating event costs (not part of the
// package). Load with Rcpp::sourceCpp("dev/robustness/weighted_events.cpp").
//
// Same greedy algorithms as greedy_cost_only() / bfb_cost_only() in
// src/fast_distances.cpp, but every block event costs
//   min(1, block_length / L)
// where block_length is the summed width (in bins) of the segments it covers.
// L = 0 (or widths all >= L) gives the usual costs.
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

static inline int sgn(double v) { return (v > 0) - (v < 0); }
static inline bool is_odd(double v) { long long r = llround(v); return std::llabs(r) % 2 == 1; }

static inline double event_cost(const std::vector<double>& cw, int head, int end, double L) {
  if (L <= 0) return 1.0;
  double len = cw[end + 1] - cw[head];          // cumulative widths
  return std::min(1.0, len / L);
}

static double greedy_w(const double* a, const double* b, int n, std::vector<double>& x,
                       const std::vector<double>& cw, double L) {
  for (int i = 0; i < n; ++i) x[i] = a[i] - b[i];
  double cost = 0; int head = 0;
  while (true) {
    while (head < n && x[head] == 0) ++head;
    if (head >= n) break;
    int s = (x[head] > 0.0) ? 1 : -1, end = head;
    while (end + 1 < n && sgn(x[end + 1]) == s) ++end;
    double step = std::abs(x[head]);
    for (int i = head + 1; i <= end; ++i) step = std::min(step, std::abs(x[i]));
    const double delta = (s < 0) ? step : -step;
    for (int i = head; i <= end; ++i) x[i] += delta;
    cost += event_cost(cw, head, end, L);
  }
  return cost;
}

static double bfb_w(const double* a, const double* b, int n, std::vector<double>& x,
                    const std::vector<double>& cw, double L) {
  for (int i = 0; i < n; ++i) x[i] = a[i] - b[i];
  double cost = 0; int head = 0;
  while (true) {
    while (head < n && !is_odd(x[head])) ++head;
    if (head >= n) break;
    int s = sgn(x[head]), end = head;
    while (end + 1 < n && sgn(x[end + 1]) == s) ++end;
    const double delta = (s < 0) ? 1.0 : -1.0;
    for (int i = head; i <= end; ++i) x[i] += delta;
    cost += event_cost(cw, head, end, L);
  }
  for (int i = 0; i < n; ) {
    if (x[i] != 0) { int j = i; while (j < n && x[j] != 0) ++j; cost += event_cost(cw, i, j - 1, L); i = j; }
    else ++i;
  }
  return cost;
}

// Pairwise min(G_w, B_w) for one cells x segments matrix; `widths` are the
// segment widths in bins. Diagonal = 0.
// [[Rcpp::export]]
NumericMatrix weighted_min_distance(NumericMatrix X, NumericVector widths, double L) {
  const int m = X.nrow(), n = X.ncol();
  std::vector<double> cw(n + 1, 0.0);
  for (int j = 0; j < n; ++j) cw[j + 1] = cw[j] + widths[j];
  std::vector<double> rows((size_t)m * n);
  for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) rows[(size_t)i * n + j] = X(i, j);
  NumericMatrix D(m, m);
  std::vector<double> buf(n);
  for (int i = 0; i < m; ++i) for (int j = i + 1; j < m; ++j) {
    const double* a = &rows[(size_t)i * n]; const double* b = &rows[(size_t)j * n];
    double d = std::min(greedy_w(a, b, n, buf, cw, L), bfb_w(a, b, n, buf, cw, L));
    D(i, j) = D(j, i) = d;
  }
  return D;
}
