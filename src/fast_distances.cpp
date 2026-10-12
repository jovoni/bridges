
#include <Rcpp.h>
using namespace Rcpp;

// helper: integer-odd check via rounding (assumes integer-like inputs)
inline bool is_odd(double v) {
  long long r = llround(v);
  return std::llabs(r) % 2 == 1;
}

inline int sgn(double v) {
  if (v > 0) return 1;
  if (v < 0) return -1;
  return 0;
}

//' Greedy block-adjustment distance (no history)
//'
//' See \code{find_greedy_distance_with_steps_fast()} for the user-facing wrapper.
//'
//' @param a NumericVector.
//' @param b NumericVector.
//' @param target_val Numeric scalar, default 2.0.
//' @return A list; see R wrapper documentation for fields.
//' @noRd
// [[Rcpp::export]]
 List greedy_distance_cpp(NumericVector a,
                          NumericVector b,
                          double target_val = 2.0) {
   const int n = a.size();
  NumericVector current = clone(a);
  NumericVector x(n);

  // x = current - b
  for (int i = 0; i < n; ++i) x[i] = current[i] - b[i];

  // Best ancestor tracking (distance to target and step index)
  double best_dist = 0.0;
  for (int i = 0; i < n; ++i) best_dist += std::abs(current[i] - target_val);
  NumericVector best_ancestor = clone(current);
  int best_step = 0;        // number of steps taken when best ancestor was seen (0 = at start)

  int cost = 0;             // total number of steps
  double cur_dist = best_dist;  // distance of `current` to target, kept up to date
  int head = 0;

  while (true) {
    // find first differing position. Positions before `head` are already
    // zero and are never modified again, so the scan resumes from there.
    while (head < n && x[head] == 0) ++head;
    if (head >= n) break; // done

    // contiguous block with same sign
    int sign_val = (x[head] > 0.0) ? 1 : -1;
    int end = head;
    while (end + 1 < n) {
      double xi = x[end + 1];
      int s = (xi > 0.0) ? 1 : ((xi < 0.0) ? -1 : 0);
      if (s == sign_val) ++end; else break;
    }

    // step size = min(abs(x[head..end]))
    double step = std::abs(x[head]);
    for (int i = head + 1; i <= end; ++i) {
      double ax = std::abs(x[i]);
      if (ax < step) step = ax;
    }

    // apply
    const double delta = (sign_val < 0) ? step : -step; // x<0 => add; x>0 => subtract
    for (int i = head; i <= end; ++i) {
      cur_dist -= std::abs(current[i] - target_val);
      current[i] += delta;
      cur_dist += std::abs(current[i] - target_val);
      x[i] = current[i] - b[i];
    }

    ++cost;

    // update best ancestor if strictly better (which.min picks first minimum)
    double d = cur_dist;
    if (d < best_dist) {
      best_dist = d;
      best_ancestor = clone(current);
      best_step = cost; // ancestor seen after this many steps from start
    }
  }

  // costs and deltas
  int a_cost = best_step;
  int b_cost = cost - best_step;

  NumericVector a_delta(n), b_delta(n);
  for (int i = 0; i < n; ++i) {
    a_delta[i] = a[i] - best_ancestor[i];
    b_delta[i] = b[i] - best_ancestor[i];
  }

  return List::create(
    _["cost"]     = cost,
    _["ancestor"] = best_ancestor,
    _["b_cost"]   = b_cost,
    _["a_cost"]   = a_cost,
    _["a_delta"]  = a_delta,
    _["b_delta"]  = b_delta
  );
 }

//' Greedy A-contiguous distance (no history)
//'
//' See \code{find_greedy_distance_A_contig_fast()} for the user-facing wrapper.
//'
//' @param a NumericVector.
//' @param b NumericVector.
//' @param penalty Numeric scalar, default 0.0.
//' @return A list; see R wrapper documentation for fields.
//' @noRd
// [[Rcpp::export]]
 List greedy_bfb_distance_cpp(NumericVector a,
                       NumericVector b,
                       double penalty = 0.0) {
   const int n = a.size();
   NumericVector current = clone(a);
   NumericVector x(n);

   // x = current - b
   for (int i = 0; i < n; ++i) x[i] = current[i] - b[i];

   double cost = 0.0;

   // loop until all differences are even (or zero). Positions before `head`
   // are even and are never modified again, so the scan resumes from there.
   int head = 0;
   while (true) {
     while (head < n && !is_odd(x[head])) ++head;
     if (head >= n) break; // all even -> done

     int sign_val = sgn(x[head]);
     int end = head;
     while (end + 1 < n && sgn(x[end + 1]) == sign_val) ++end;

     // +/-1 to entire block
     const double delta = (sign_val < 0) ? 1.0 : -1.0;
     for (int i = head; i <= end; ++i) {
       current[i] += delta;
       x[i] = current[i] - b[i];
     }
     cost += 1.0;
   }

   // count number of TRUE runs in rle(x != 0)
   int nonzero_runs = 0;
   int i = 0;
   while (i < n) {
     if (x[i] != 0.0) {
       ++nonzero_runs;
       // advance through the whole nonzero run
       int j = i + 1;
       while (j < n && x[j] != 0.0) ++j;
       i = j;
     } else {
       ++i;
     }
   }

   cost += static_cast<double>(nonzero_runs) + penalty;

   // ancestor is average of final state and target
   NumericVector ancestor(n), a_delta(n), b_delta(n);
   for (int k = 0; k < n; ++k) {
     ancestor[k] = 0.5 * (current[k] + b[k]);
     a_delta[k] = a[k] - ancestor[k];
     b_delta[k] = b[k] - ancestor[k];
   }

   // costs as in your R version
   double a_cost = cost - 1.0;
   double b_cost = 1.0;

   return List::create(
     _["history"]  = R_NilValue,      // no history stored
     _["cost"]     = cost,
     _["ancestor"] = ancestor,
     _["b_cost"]   = b_cost,
     _["a_cost"]   = a_cost,
     _["a_delta"]  = a_delta,
     _["b_delta"]  = b_delta,
     _["left"]     = current,
     _["right"]    = b
   );
 }

// ── Whole-matrix distances ───────────────────────────────────────────────────
//
// Cost-only versions of greedy_distance_cpp() and greedy_bfb_distance_cpp():
// same algorithms, no ancestor or delta vectors. `x` is a work buffer of
// length n that receives a - b.

static double greedy_cost_only(const double* a, const double* b, int n,
                               std::vector<double>& x) {
  for (int i = 0; i < n; ++i) x[i] = a[i] - b[i];
  int cost = 0;
  int head = 0;
  while (true) {
    while (head < n && x[head] == 0) ++head;
    if (head >= n) break;
    int sign_val = (x[head] > 0.0) ? 1 : -1;
    int end = head;
    while (end + 1 < n && sgn(x[end + 1]) == sign_val) ++end;
    double step = std::abs(x[head]);
    for (int i = head + 1; i <= end; ++i) step = std::min(step, std::abs(x[i]));
    const double delta = (sign_val < 0) ? step : -step;
    for (int i = head; i <= end; ++i) x[i] += delta;
    ++cost;
  }
  return cost;
}

static double bfb_cost_only(const double* a, const double* b, int n,
                            std::vector<double>& x, double penalty) {
  for (int i = 0; i < n; ++i) x[i] = a[i] - b[i];
  double cost = 0.0;
  int head = 0;
  while (true) {
    while (head < n && !is_odd(x[head])) ++head;
    if (head >= n) break;
    int sign_val = sgn(x[head]);
    int end = head;
    while (end + 1 < n && sgn(x[end + 1]) == sign_val) ++end;
    const double delta = (sign_val < 0) ? 1.0 : -1.0;
    for (int i = head; i <= end; ++i) x[i] += delta;
    cost += 1.0;
  }
  int nonzero_runs = 0;
  for (int i = 0; i < n; ) {
    if (x[i] != 0.0) {
      ++nonzero_runs;
      while (i < n && x[i] != 0.0) ++i;
    } else {
      ++i;
    }
  }
  return cost + nonzero_runs + penalty;
}

//' Pairwise greedy and BFB distance matrices for one cells x segments matrix
//'
//' Computes, for every pair of rows of \code{X}, the cost of
//' \code{greedy_distance_cpp()} (G) and, if \code{use_bfb}, of
//' \code{greedy_bfb_distance_cpp()} (B), and their element-wise minimum.
//' Identical rows are detected first and distances are computed once per
//' pair of distinct rows. Diagonals are set to \code{Inf}, as in
//' \code{compute_distance_matrix()}.
//'
//' @param X Numeric matrix, cells x segments.
//' @param use_bfb Whether to compute the BFB distance; if \code{FALSE},
//'   \code{min} equals the greedy distance.
//' @param penalty Penalty added to every BFB cost.
//' @param keep_parts Whether to also return the full G and B matrices.
//' @return A list with \code{min}, \code{greedy} and \code{bfb} (the last two
//'   \code{NULL} unless \code{keep_parts}), and \code{n_unique}, the number of
//'   distinct rows.
//' @noRd
// [[Rcpp::export]]
List distance_matrices_cpp(NumericMatrix X, bool use_bfb, double penalty = 0.0,
                           bool keep_parts = false) {
  const int m = X.nrow(), n = X.ncol();

  // Row-major copy so each cell's profile is contiguous
  std::vector<double> rows((size_t)m * n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j) rows[(size_t)i * n + j] = X(i, j);
  auto row = [&](int i) { return rows.data() + (size_t)i * n; };

  // Group identical rows: sort row indices lexicographically, then assign
  // one id per run of equal rows.
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

  // Distances between distinct rows (u x u, row-major). Identical rows have
  // G = 0 and B = penalty, exactly what the kernels return for a == b.
  std::vector<double> G((size_t)u * u, 0.0), B;
  if (use_bfb) B.assign((size_t)u * u, penalty);
  std::vector<double> buf(n);
  for (int p = 0; p < u; ++p) {
    for (int q = p + 1; q < u; ++q) {
      const double* a = row(rep[p]);
      const double* b = row(rep[q]);
      double g = greedy_cost_only(a, b, n, buf);
      G[(size_t)p * u + q] = G[(size_t)q * u + p] = g;
      if (use_bfb) {
        double bb = bfb_cost_only(a, b, n, buf, penalty);
        B[(size_t)p * u + q] = B[(size_t)q * u + p] = bb;
      }
    }
  }

  // Expand to all rows
  NumericMatrix Dmin(m, m);
  NumericMatrix Dg = keep_parts ? NumericMatrix(m, m) : NumericMatrix(0, 0);
  NumericMatrix Db = (keep_parts && use_bfb) ? NumericMatrix(m, m) : NumericMatrix(0, 0);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      size_t k = (size_t)group[i] * u + group[j];
      double g = G[k];
      double bb = use_bfb ? B[k] : R_PosInf;
      if (i == j) { g = R_PosInf; bb = R_PosInf; }
      Dmin(i, j) = std::min(g, bb);
      if (keep_parts) {
        Dg(i, j) = g;
        if (use_bfb) Db(i, j) = bb;
      }
    }
  }

  return List::create(
    _["min"]      = Dmin,
    _["greedy"]   = keep_parts ? (SEXP)Dg : R_NilValue,
    _["bfb"]      = (keep_parts && use_bfb) ? (SEXP)Db : R_NilValue,
    _["n_unique"] = u
  );
}
