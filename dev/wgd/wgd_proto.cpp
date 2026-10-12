// Prototype WGD-aware distances (not part of the package).
// Load with Rcpp::sourceCpp("dev/wgd/wgd_proto.cpp").
//
// W(a -> b, f): cost of explaining b as f * x plus greedy events, where x is
// reached from a by greedy events:  min_x G(a, x) + G(f * x, b).
// The candidates for x are all intermediate profiles on the greedy path from
// a towards floor(b / f) and towards ceil(b / f) (both endpoints included).
// This is an upper bound on the exact minimum.
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

static inline int sgn(double v) { return (v > 0) - (v < 0); }
static inline bool is_odd(double v) { long long r = llround(v); return std::llabs(r) % 2 == 1; }

// Same as greedy_cost_only() in src/fast_distances.cpp
static double greedy_cost(const double* a, const double* b, int n, std::vector<double>& x) {
  for (int i = 0; i < n; ++i) x[i] = a[i] - b[i];
  int cost = 0, head = 0;
  while (true) {
    while (head < n && x[head] == 0) ++head;
    if (head >= n) break;
    int s = (x[head] > 0.0) ? 1 : -1, end = head;
    while (end + 1 < n && sgn(x[end + 1]) == s) ++end;
    double step = std::abs(x[head]);
    for (int i = head + 1; i <= end; ++i) step = std::min(step, std::abs(x[i]));
    const double delta = (s < 0) ? step : -step;
    for (int i = head; i <= end; ++i) x[i] += delta;
    ++cost;
  }
  return cost;
}

// Same as bfb_cost_only() in src/fast_distances.cpp
static double bfb_cost(const double* a, const double* b, int n, std::vector<double>& x) {
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

// Walk the greedy path from a to target; at every state x (including a and
// target) evaluate steps + G(f * x, b). Return the minimum.
static double wgd_path_cost(const double* a, const double* target, const double* b,
                            int n, double f, std::vector<double>& cur,
                            std::vector<double>& diff, std::vector<double>& fx,
                            std::vector<double>& buf,
                            std::vector<double>* best_x = nullptr, int* best_steps = nullptr) {
  for (int i = 0; i < n; ++i) { cur[i] = a[i]; diff[i] = a[i] - target[i]; }
  auto eval = [&](int steps) {
    for (int i = 0; i < n; ++i) fx[i] = f * cur[i];
    return steps + greedy_cost(fx.data(), b, n, buf);
  };
  int steps = 0, head = 0;
  double best = eval(0);
  if (best_x) { *best_x = cur; *best_steps = 0; }
  while (true) {
    while (head < n && diff[head] == 0) ++head;
    if (head >= n) break;
    int s = (diff[head] > 0.0) ? 1 : -1, end = head;
    while (end + 1 < n && sgn(diff[end + 1]) == s) ++end;
    double step = std::abs(diff[head]);
    for (int i = head + 1; i <= end; ++i) step = std::min(step, std::abs(diff[i]));
    const double delta = (s < 0) ? step : -step;
    for (int i = head; i <= end; ++i) { cur[i] += delta; diff[i] += delta; }
    ++steps;
    if (steps >= best) break;  // remaining states cannot beat the best
    double c = eval(steps);
    if (c < best) {
      best = c;
      if (best_x) { *best_x = cur; *best_steps = steps; }
    }
  }
  return best;
}

static double wgd_cost_impl(const double* a, const double* b, int n, double f,
                            std::vector<double>& t, std::vector<double>& cur,
                            std::vector<double>& diff, std::vector<double>& fx,
                            std::vector<double>& buf) {
  for (int i = 0; i < n; ++i) t[i] = std::floor(b[i] / f);
  double best = wgd_path_cost(a, t.data(), b, n, f, cur, diff, fx, buf);
  bool same = true;
  for (int i = 0; i < n; ++i) { double c = std::ceil(b[i] / f); if (c != t[i]) same = false; t[i] = c; }
  if (!same) best = std::min(best, wgd_path_cost(a, t.data(), b, n, f, cur, diff, fx, buf));
  return best;
}

// Best ancestor for W(a -> b, f): returns cost, x, steps (events a -> x) and
// events after the doubling (f * x -> b).
// [[Rcpp::export]]
List wgd_ancestor(NumericVector a, NumericVector b, double f) {
  int n = a.size();
  std::vector<double> t(n), cur(n), diff(n), fx(n), buf(n), x1, x2;
  int s1 = 0, s2 = 0;
  for (int i = 0; i < n; ++i) t[i] = std::floor(b[i] / f);
  double c1 = wgd_path_cost(a.begin(), t.data(), b.begin(), n, f, cur, diff, fx, buf, &x1, &s1);
  for (int i = 0; i < n; ++i) t[i] = std::ceil(b[i] / f);
  double c2 = wgd_path_cost(a.begin(), t.data(), b.begin(), n, f, cur, diff, fx, buf, &x2, &s2);
  bool first = c1 <= c2;
  double cost = first ? c1 : c2;
  return List::create(_["cost"] = cost,
                      _["x"] = NumericVector((first ? x1 : x2).begin(), (first ? x1 : x2).end()),
                      _["steps"] = first ? s1 : s2,
                      _["after"] = cost - (first ? s1 : s2));
}

// [[Rcpp::export]]
double wgd_cost(NumericVector a, NumericVector b, double f) {
  int n = a.size();
  std::vector<double> t(n), cur(n), diff(n), fx(n), buf(n);
  return wgd_cost_impl(a.begin(), b.begin(), n, f, t, cur, diff, fx, buf);
}

// For one cells x segments matrix: S0 = pmin(G, B), and directional WGD costs
// W1[i, j] = W(row i -> row j, 2), W2[i, j] = W(row i -> row j, 4).
// Diagonals are Inf. Identical rows are computed once.
// [[Rcpp::export]]
List wgd_matrices(NumericMatrix X) {
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
  std::vector<double> t(n), cur(n), diff(n), fx(n), buf(n);
  for (int p = 0; p < u; ++p) {
    for (int q = 0; q < u; ++q) {
      const double* a = row(rep[p]);
      const double* b = row(rep[q]);
      if (q > p) {
        double s0 = std::min(greedy_cost(a, b, n, buf), bfb_cost(a, b, n, buf));
        S0[(size_t)p * u + q] = S0[(size_t)q * u + p] = s0;
      }
      W1[(size_t)p * u + q] = wgd_cost_impl(a, b, n, 2.0, t, cur, diff, fx, buf);
      W2[(size_t)p * u + q] = wgd_cost_impl(a, b, n, 4.0, t, cur, diff, fx, buf);
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
