
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

  while (true) {
    // find first differing position
    int head = -1;
    for (int i = 0; i < n; ++i) {
      if (x[i] != 0) { head = i; break; }
    }
    if (head < 0) break; // done

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
      current[i] += delta;
      x[i] = current[i] - b[i];
    }

    ++cost;

    // update best ancestor if strictly better (which.min picks first minimum)
    double d = 0.0;
    for (int i = 0; i < n; ++i) d += std::abs(current[i] - target_val);
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

   // loop until all differences are even (or zero)
   while (true) {
     int head = -1;
     for (int i = 0; i < n; ++i) {
       if (is_odd(x[i])) { head = i; break; }
     }
     if (head < 0) break; // all even -> done

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
