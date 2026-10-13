// Prototype event-model distances (not part of the package).
// Load with Rcpp::sourceCpp("dev/wgd/event_model.cpp").
//
// Event model: every event changes the copy number of one contiguous segment
// by exactly +1 or -1 (no restriction at copy number 0). WGD multiplies the
// whole profile by 2.
//
// Without WGD, the minimum number of events between a and b is
//   E(a, b) = TV(b - a) / 2,   TV(v) = |v_1| + sum_i |v_i - v_{i-1}| + |v_n|
// (each event changes TV by at most 2, and any pair of opposite jumps can be
// removed by one event).
//
// General ancestor problem, solved exactly by dynamic programming over bins
// (the state is the integer value of x in each bin):
//   min_x  E(fp * x, p) + E(fq * x, q)
// - fp = fq = 1: ancestor of two cells at the same WGD level
// - fp = 1, fq = 2 (or 4): b descends from a doubling of x
// - fp = fq = 2: both children doubled after the ancestor
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

static double tv_diff(const double* p, double fp, const double* x, int n) {
  // TV(p - fp * x) with zero padding
  double tv = 0.0, prev = 0.0;
  for (int i = 0; i < n; ++i) {
    double v = p[i] - fp * x[i];
    tv += std::abs(v - prev);
    prev = v;
  }
  return tv + std::abs(prev);
}

static double event_cost(const double* a, const double* b, int n) {
  double tv = 0.0, prev = 0.0;
  for (int i = 0; i < n; ++i) {
    double v = b[i] - a[i];
    tv += std::abs(v - prev);
    prev = v;
  }
  return (tv + std::abs(prev)) / 2.0;
}

// Exact minimiser of TV(p - fp x) + TV(q - fq x), tie-broken towards x = t.
// Writes the best x to x_out and returns the event cost (TV sum / 2).
static double ancestor_dp(const double* p, double fp, const double* q, double fq, int n,
                          double t, std::vector<double>& x_out,
                          std::vector<double>& cost, std::vector<double>& next,
                          std::vector<int>& back) {
  double maxv = 0.0;
  for (int i = 0; i < n; ++i) maxv = std::max(maxv, std::max(p[i] / fp, q[i] / fq));
  const int K = (int)std::ceil(maxv) + 2;            // x in 0 .. K-1
  const double gamma = 1.0 / (2.0 * n * (K + std::abs(t)) + 2.0);  // tie-break only
  cost.assign(K, 0.0);
  next.assign(K, 0.0);
  back.assign((size_t)n * K, 0);

  for (int x = 0; x < K; ++x)
    cost[x] = std::abs(p[0] - fp * x) + std::abs(q[0] - fq * x) + gamma * std::abs(x - t);
  for (int i = 1; i < n; ++i) {
    for (int x = 0; x < K; ++x) {
      double best = R_PosInf; int arg = 0;
      double vp = p[i] - fp * x, vq = q[i] - fq * x;
      for (int y = 0; y < K; ++y) {
        double c = cost[y] + std::abs(vp - (p[i - 1] - fp * y)) + std::abs(vq - (q[i - 1] - fq * y));
        if (c < best) { best = c; arg = y; }
      }
      next[x] = best + gamma * std::abs(x - t);
      back[(size_t)i * K + x] = arg;
    }
    std::swap(cost, next);
  }
  double best = R_PosInf; int arg = 0;
  for (int x = 0; x < K; ++x) {
    double c = cost[x] + std::abs(p[n - 1] - fp * x) + std::abs(q[n - 1] - fq * x);
    if (c < best) { best = c; arg = x; }
  }
  x_out.assign(n, 0.0);
  for (int i = n - 1; i >= 0; --i) {
    x_out[i] = arg;
    if (i > 0) arg = back[(size_t)i * K + arg];
  }
  // exact event cost of the chosen x (without the tie-break term)
  return (tv_diff(p, fp, x_out.data(), n) + tv_diff(q, fq, x_out.data(), n)) / 2.0;
}

static inline int sgn(double v) { return (v > 0) - (v < 0); }
static inline bool is_odd(double v) { long long r = llround(v); return std::llabs(r) % 2 == 1; }

// BFB distance, same algorithm as bfb_cost_only() in src/fast_distances.cpp
static double bfb_cost(const double* a, const double* b, int n, std::vector<double>& x) {
  x.resize(n);
  for (int i = 0; i < n; ++i) x[i] = a[i] - b[i];
  double cost = 0; int head = 0;
  while (true) {
    while (head < n && !is_odd(x[head])) ++head;
    if (head >= n) break;
    int s = sgn(x[head]), end = head;
    while (end + 1 < n && sgn(x[end + 1]) == s) ++end;
    const double delta = (s < 0) ? 1.0 : -1.0;
    for (int i = head; i <= end; ++i) x[i] += delta;
    cost += 1.0;
  }
  int runs = 0;
  for (int i = 0; i < n; ) { if (x[i] != 0) { ++runs; while (i < n && x[i] != 0) ++i; } else ++i; }
  return cost + runs;
}

// [[Rcpp::export]]
double event_distance(NumericVector a, NumericVector b) {
  return event_cost(a.begin(), b.begin(), a.size());
}

// Ancestor x minimising E(fp x, p) + E(fq x, q); tie-break towards target t.
// Returns cost, x, and the two branch costs.
// [[Rcpp::export]]
List event_ancestor(NumericVector p, double fp, NumericVector q, double fq, double t = 1.0) {
  int n = p.size();
  std::vector<double> x, c1, c2; std::vector<int> back;
  double cost = ancestor_dp(p.begin(), fp, q.begin(), fq, n, t, x, c1, c2, back);
  return List::create(_["cost"] = cost,
                      _["x"] = NumericVector(x.begin(), x.end()),
                      _["cost_p"] = tv_diff(p.begin(), fp, x.data(), n) / 2.0,
                      _["cost_q"] = tv_diff(q.begin(), fq, x.data(), n) / 2.0);
}

// Same layout as wgd_matrices() in wgd_proto.cpp, under the event model:
// S0 = E(a, b) (or min(E, B)); W1[i, j] = min_x E(x, a_i) + E(2x, a_j); W2 the same with 4x.
// Diagonals are Inf. Identical rows are computed once.
// With use_bfb, S0 = min(E, B) per chromosome/allele, as fit() does with G and B.
// [[Rcpp::export]]
List event_matrices(NumericMatrix X, bool use_bfb = false) {
  const int m = X.nrow(), n = X.ncol();
  std::vector<double> rows((size_t)m * n);
  for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) rows[(size_t)i * n + j] = X(i, j);
  auto row = [&](int i) { return rows.data() + (size_t)i * n; };

  std::vector<int> order(m);
  for (int i = 0; i < m; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](int p, int q) {
    return std::lexicographical_compare(row(p), row(p) + n, row(q), row(q) + n);
  });
  std::vector<int> group(m), rep;
  for (int k = 0; k < m; ++k) {
    int i = order[k];
    if (k == 0 || !std::equal(row(i), row(i) + n, row(order[k - 1]))) rep.push_back(i);
    group[i] = (int)rep.size() - 1;
  }
  const int u = (int)rep.size();

  std::vector<double> S0((size_t)u * u, 0.0), W1((size_t)u * u), W2((size_t)u * u);
  std::vector<double> x, c1, c2; std::vector<int> back;
  for (int p = 0; p < u; ++p) {
    for (int q = 0; q < u; ++q) {
      const double* a = row(rep[p]);
      const double* b = row(rep[q]);
      if (q > p) {
        double s0 = event_cost(a, b, n);
        if (use_bfb) s0 = std::min(s0, bfb_cost(a, b, n, c1));
        S0[(size_t)p * u + q] = S0[(size_t)q * u + p] = s0;
      }
      W1[(size_t)p * u + q] = ancestor_dp(a, 1.0, b, 2.0, n, 0.0, x, c1, c2, back);
      W2[(size_t)p * u + q] = ancestor_dp(a, 1.0, b, 4.0, n, 0.0, x, c1, c2, back);
    }
  }

  NumericMatrix M0(m, m), M1(m, m), M2(m, m);
  for (int j = 0; j < m; ++j) for (int i = 0; i < m; ++i) {
    size_t k = (size_t)group[i] * u + group[j];
    bool d = (i == j);
    M0(i, j) = d ? R_PosInf : S0[k];
    M1(i, j) = d ? R_PosInf : W1[k];
    M2(i, j) = d ? R_PosInf : W2[k];
  }
  return List::create(_["S0"] = M0, _["W1"] = M1, _["W2"] = M2, _["n_unique"] = u);
}
